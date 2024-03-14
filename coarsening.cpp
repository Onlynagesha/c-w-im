#include "coarsening.h"
#include "dump.h"
#include "utils/easylog.h"
#include "wim.h"
#include <fmt/ranges.h>
#include <nwgraph/adaptors/edge_range.hpp>

namespace {
constexpr auto MAX_N_MEMBERS = CoarsenedVertexDetails::MAX_N_MEMBERS;
constexpr auto P_NOT_ASSIGNED = -1.0_ep;

using VertexPair = std::pair<vertex_id_t, vertex_id_t>;

auto merge_p(edge_probability_t e1, edge_probability_t e2, NeighborMatchRule rule) -> edge_probability_t {
  switch (rule) {
  case NeighborMatchRule::HEM_P_MAX:
  case NeighborMatchRule::LEM_P_MAX:
    return std::max(e1, e2);
  case NeighborMatchRule::HEM_P_PRODUCT:
  case NeighborMatchRule::LEM_P_PRODUCT:
    return 1.0_ep - (1.0_ep - e1) * (1.0_ep - e2);
  default:
    BOOST_ASSERT(!"Unreachable branch.");
    std::abort();
  }
}

template <is_edge_property E>
auto merge_p(const E& e1, const E& e2, NeighborMatchRule rule) -> edge_probability_t {
  return merge_p(e1.p, e2.p, rule);
}

auto p_infinity_by_rule(NeighborMatchRule rule) -> edge_probability_t {
  switch (rule) {
  case NeighborMatchRule::HEM_P_MAX:
  case NeighborMatchRule::HEM_P_PRODUCT:
    return std::numeric_limits<edge_probability_t>::lowest(); // -inf
  case NeighborMatchRule::LEM_P_MAX:
  case NeighborMatchRule::LEM_P_PRODUCT:
    return std::numeric_limits<edge_probability_t>::max(); // +inf
  default:
    BOOST_ASSERT(!"Unreachable branch.");
    std::abort();
  }
}

auto is_better_p(edge_probability_t before, edge_probability_t after, NeighborMatchRule rule) {
  switch (rule) {
  case NeighborMatchRule::HEM_P_MAX:
  case NeighborMatchRule::HEM_P_PRODUCT:
    return after > before; // HEM: higher is better
  case NeighborMatchRule::LEM_P_MAX:
  case NeighborMatchRule::LEM_P_PRODUCT:
    return after < before; // LEM: lower is better
  default:
    BOOST_ASSERT(!"Unreachable branch.");
    std::abort();
  }
}

template <size_t N, size_t... Is>
  requires(N >= 2)
constexpr auto get_next_j_helper(size_t j0, std::index_sequence<Is...>) {
  if constexpr (sizeof...(Is) == 1) {
    return (j0 + (Is, ...)) % N;
  } else {
    return std::tuple{((j0 + Is) % N)...};
  }
}

template <size_t N>
  requires(N >= 2)
constexpr auto get_next_j(size_t j0) {
  BOOST_ASSERT_MSG(j0 < N, "j0 is out of range [0, N).");
  return get_next_j_helper<N>(j0, make_index_sequence_by_offset<N - 1, 1>{});
}

template <class T>
struct VertexPairWithValue {
  vertex_id_t u;
  vertex_id_t v;
  [[no_unique_address]] T value;

  auto vertices() const {
    return std::tuple{u, v};
  }
};

template <class T>
auto radix_sort_vertex_pairs(std::span<VertexPairWithValue<T>> vertex_pairs, vertex_id_t n) -> void {
  auto m = vertex_pairs.size();
  auto count = std::vector<vertex_id_t>(n);
  auto temp = std::vector<VertexPairWithValue<T>>(m);
  MYLOG_FMT_TRACE("vertex_pairs before radix sorting: {}", vertex_pairs | TRANSFORM_VIEW(std::pair{_1.u, _1.v}));
  // Step 1: Sorts by v, vertex_pairs -> temp
  for (const auto& item : vertex_pairs) {
    BOOST_ASSERT_MSG(item.u < n && item.v < n, "v is out of range [0, n).");
    count[item.v] += 1;
  }
  std::exclusive_scan(count.begin(), count.end(), count.begin(), 0);
  for (const auto& item : vertex_pairs) {
    temp[count[item.v]] = item;
    count[item.v] += 1;
  }
  MYLOG_FMT_TRACE("vertex_pairs after radix sorting step 1: {}", temp | TRANSFORM_VIEW(std::pair{_1.u, _1.v}));
  // Step 2: Sorts by u, temp -> vertex_pairs
  ranges::fill(count, 0);
  for (const auto& item : temp) {
    count[item.u] += 1;
  }
  std::exclusive_scan(count.begin(), count.end(), count.begin(), 0);
  for (const auto& item : temp) {
    vertex_pairs[count[item.u]] = item;
    count[item.u] += 1;
  }
  MYLOG_FMT_TRACE("vertex_pairs after radix sorting step 2: {}", //
                  vertex_pairs | TRANSFORM_VIEW(std::pair{_1.u, _1.v}));
  BOOST_ASSERT_MSG(ranges::is_sorted(vertex_pairs | TRANSFORM_VIEW(std::pair{_1.u, _1.v})), //
                   "Implementation error: Not sorted correctly.");
}

template <class Range>
auto radix_sort_vertex_pairs(Range& vertex_pairs, vertex_id_t n) {
  radix_sort_vertex_pairs(std::span{vertex_pairs.begin(), vertex_pairs.end()}, n);
}

template <is_edge_property E>
auto merge_edge_to_undirected_generic(AdjacencyList<E>& graph, const CoarseningParams& params)
    -> AdjacencyList<edge_probability_t> {
  using ItemType = VertexPairWithValue<edge_probability_t>;
  auto e_pairs = make_reserved_vector<ItemType>(graph.num_edges() + 1);
  for (auto [u, v, w] : graph::make_edge_range<0>(graph)) {
    if (u == v) {
      continue; // Skips self-loops
    }
    (u < v) ? e_pairs.push_back({u, v, w.p}) : e_pairs.push_back({v, u, w.p});
  }
  radix_sort_vertex_pairs(e_pairs, graph::num_vertices(graph));
  e_pairs.push_back({static_cast<vertex_id_t>(-1), static_cast<vertex_id_t>(-1), 0.0_ep}); // Dummy end

  constexpr auto check_duplication = [](std::span<const ItemType> items) {
    auto chunks = items | views::chunk_by(LAMBDA_2(_1.u == _2.u && _1.v == _2.v));
    return ranges::max(chunks | views::transform(ranges::size)) <= 2;
  };
  BOOST_ASSERT_MSG(check_duplication(e_pairs), "Unexpected input data: duplicated edge detected.");

  // Ensures |V| is unchanged
  auto edge_list = UndirectedEdgeList<edge_probability_t>{graph::num_vertices(graph)};
  edge_list.open_for_push_back();
  auto e_pairs_adj = e_pairs | views::adjacent<2>;
  for (auto it = e_pairs_adj.begin(); it != e_pairs_adj.end();) {
    const auto& [p1, p2] = *it;
    if (p1.u == p2.u && p1.v == p2.v) {
      edge_list.push_back(p1.u, p1.v, merge_p(p1.value, p2.value, params.neighbor_match_rule));
      it += 2; // Merges (p1, p2)
    } else {
      edge_list.push_back(p1.u, p1.v, p1.value);
      it += 1; // p1 only
    }
  }
  edge_list.close_for_push_back();
  return AdjacencyList<edge_probability_t>{edge_list};
}

// vertex_weights can be either a random-access container or something dummy like views::repeat
template <auto WhichP,                   // E::*WhichP -> edge_probability_t
          is_edge_property E, int IsInv, // IsInv == 1 is transposed graph; IsInv == 0 otherwise
          random_access_range_of<vertex_weight_t> VertexWeights>
  requires(member_object_pointer_from_to<decltype(WhichP), E, edge_probability_t> && (IsInv == 0 || IsInv == 1))
auto get_heuristic_generic(const graph::adjacency<IsInv, E>& graph, std::span<const vertex_id_t> group_id,
                           VertexWeights&& vertex_weights, vertex_id_t v, InOutHeuristicRule rule)
    -> edge_probability_t {
  auto res = 0.0_ep;
  for (auto [u, w] : graph[v]) {
    if (group_id[u] == group_id[v]) {
      continue;
    }
    if (rule == InOutHeuristicRule::P) {
      res += w.*WhichP;
    } else if (rule == InOutHeuristicRule::W) {
      res += w.*WhichP * vertex_weights[u];
    } else {
      res += 1.0_ep; // Counts only
    }
  }
  return res;
}

template <is_edge_property E, random_access_range_of<vertex_weight_t> VertexWeights>
auto get_heuristic_in(const InvAdjacencyList<E>& inv_graph, std::span<const vertex_id_t> group_id,
                      VertexWeights&& vertex_weights, vertex_id_t v, InOutHeuristicRule rule) -> edge_probability_t {
  return get_heuristic_generic<&E::p>(inv_graph, group_id, std::forward<VertexWeights>(vertex_weights), v, rule);
}

template <is_edge_property E, random_access_range_of<vertex_weight_t> VertexWeights>
auto get_heuristic_out(const AdjacencyList<E>& graph, std::span<const vertex_id_t> group_id,
                       VertexWeights&& vertex_weights, vertex_id_t v, InOutHeuristicRule rule) -> edge_probability_t {
  return get_heuristic_generic<&E::p>(graph, group_id, std::forward<VertexWeights>(vertex_weights), v, rule);
}

template <is_edge_property E, random_access_range_of<vertex_weight_t> VertexWeights>
auto get_heuristic_out_seed(const AdjacencyList<E>& graph, std::span<const vertex_id_t> group_id,
                            VertexWeights&& vertex_weights, vertex_id_t v, InOutHeuristicRule rule)
    -> edge_probability_t {
  return get_heuristic_generic<&E::p_seed>(graph, group_id, std::forward<VertexWeights>(vertex_weights), v, rule);
}

struct GetHeuristicPResult {
  CoarsenedVertexDetails::HeuristicsContainer p;
  bool is_all_zero;

  auto normalize(InOutHeuristicRule rule) -> void {
    constexpr auto EPS = 1.0e-8_ep;
    // Ensures that each value falls in range [0, +inf) in case of floating point error
    for (auto& h_item : p) {
      BOOST_ASSERT(!std::isinf(h_item) && !std::isnan(h_item));
      h_item = std::max(h_item, 0.0_ep);
    }
    // For UNIT rule, normalizes to {0, 1} first.
    if (rule == InOutHeuristicRule::UNIT) {
      for (auto& h_item : p) {
        h_item = (h_item < EPS) ? 0.0_ep : 1.0_ep;
      }
    }
    // Normalizes by sum such that accumulate_sum(p_in) == 1.0
    auto sum = accumulate_sum(p);
    if (sum <= 0.0) {
      is_all_zero = true;
    } else {
      ranges::for_each(p, LAMBDA_1(_1 /= sum));
      is_all_zero = false;
    }
  }
};

auto get_normalized_heuristic_p_in(const CoarsenedVertexDetails& gu, InOutHeuristicRule rule) -> GetHeuristicPResult {
  auto res = GetHeuristicPResult{.p = gu.heuristics_in, .is_all_zero = false};
  res.normalize(rule);
  return res;
}

auto get_normalized_heuristic_p_out(const CoarsenedVertexDetails& gu, InOutHeuristicRule rule) -> GetHeuristicPResult {
  auto res = GetHeuristicPResult{.p = gu.heuristics_out, .is_all_zero = false};
  res.normalize(rule);
  return res;
}

auto get_normalized_heuristic_p_out_seed(const CoarsenedVertexDetails& gu, InOutHeuristicRule rule)
    -> GetHeuristicPResult {
  auto res = GetHeuristicPResult{.p = gu.heuristics_out_seed, .is_all_zero = false};
  res.normalize(rule);
  return res;
}

auto get_normalized_heuristic_p_in(const CoarsenedEdgeDetails& edge, const CoarsenedEdgeDetails* inv_edge,
                                   const CoarsenedVertexDetails& gu, const CoarsenedVertexDetails& gv,
                                   InOutHeuristicRule rule) -> GetHeuristicPResult {
  constexpr auto EPS = 1.0e-8_ep;
  auto res = GetHeuristicPResult{.p = gu.heuristics_in, .is_all_zero = false};
  if (inv_edge != nullptr) {
    // Removes the directed edge gu.members[j] <-- gv.members[i]
    for (auto [j, i] : views::cartesian_product(gu.member_indices(), gv.member_indices())) {
      if (rule == InOutHeuristicRule::P) {
        res.p[j] -= inv_edge->p_cross[i][j];
      } else if (rule == InOutHeuristicRule::W) {
        res.p[j] -= inv_edge->p_cross[i][j] * gv.vertex_weights[i];
      } else {
        BOOST_ASSERT(rule == InOutHeuristicRule::UNIT || rule == InOutHeuristicRule::COUNT);
        res.p[j] -= (inv_edge->p_cross[i][j] < EPS) ? 0.0_ep : 1.0_ep;
      }
    }
    MYLOG_FMT_TRACE("In-heuristics of group {} after removing back edges from {} = {}", //
                    gu.members, gv.members, res.p);
  }
  res.normalize(rule);
  return res;
}

// The probability of spread success from gu.members[j0] to any member of gv
// is_seed: Whether the source vertex gu.members[j0] is a seed.
auto get_merged_p_simple(const CoarsenedEdgeDetails& edge, const CoarsenedVertexDetails& gu,
                         const CoarsenedVertexDetails& gv, size_t j0, bool is_seed) -> edge_probability_t {
  auto [p_cross_first, p_internal_first] =
      is_seed ? std::tie(edge.p_seed_cross, gu.p_seed_internal) : std::tie(edge.p_cross, gu.p_internal);
  auto N = gu.n_members();

  if (N == 2) {
    // 2 : For each destination index i, at most 2 routes (see below)
    constexpr auto MAX_N_PATHS = 2 * MAX_N_MEMBERS;
    auto paths = StaticVector<edge_probability_t, MAX_N_PATHS>{};
    auto j1 = get_next_j<2>(j0);
    for (auto i : gv.member_indices()) {
      paths.push_back(p_cross_first[j0][i]);                           // (1) j0 -> i
      paths.push_back(p_internal_first[j0][j1] * edge.p_cross[j1][i]); // (2) j0 -> j1 -> i
    }
    return at_least_1_probability_r(paths);
  } else {
    BOOST_ASSERT_MSG(N == 3, "Group size large than 3 is not supported.");
    // 5 : For each destination index i, at most 5 routes (see below)
    constexpr auto MAX_N_PATHS = 5 * MAX_N_MEMBERS;
    auto paths = StaticVector<edge_probability_t, MAX_N_PATHS>{};
    auto [j1, j2] = get_next_j<3>(j0);
    for (auto i : gv.member_indices()) {
      paths.push_back(p_cross_first[j0][i]);                           // (1) j0 -> i
      paths.push_back(p_internal_first[j0][j1] * edge.p_cross[j1][i]); // (2) j0 -> j1 -> i
      paths.push_back(p_internal_first[j0][j2] * edge.p_cross[j2][i]); // (3) j0 -> j2 -> i
      // (4,5) j0 -> j1|j2 -> j2|j1 -> i
      paths.push_back(p_internal_first[j0][j1] * gu.p_internal[j1][j2] * edge.p_cross[j2][i]);
      paths.push_back(p_internal_first[j0][j2] * gu.p_internal[j2][j1] * edge.p_cross[j1][i]);
    }
    return at_least_1_probability_r(paths);
  }
}

// F[2^N][2^N]
using PreciseStateContainer = Array2D<edge_probability_t, (1 << MAX_N_MEMBERS), (1 << MAX_N_MEMBERS)>;

auto make_initial_precise_state_container() -> PreciseStateContainer {
  auto res = PreciseStateContainer{};
  ranges::fill(res | views::join, P_NOT_ASSIGNED);
  return res;
}

// T = Subset of gu.members that have received message in all the previous turns
// S = Subset of gu.members that receives message in the last turn.
// Both T, S are encoded as bitset. S is always a non-empty subset of T.
auto get_merged_p_precise_impl(const CoarsenedEdgeDetails& edge, const CoarsenedVertexDetails& gu,
                               const CoarsenedVertexDetails& gv, PreciseStateContainer& F, size_t T, size_t S)
    -> edge_probability_t {
  auto N = gu.n_members();
  auto U = (1zu << N) - 1; // U = Universal set
  BOOST_ASSERT_MSG(S <= U && T <= U, "S, T must be in the range [0, 2^N).");
  // Dynamic programming implemented as memorized searching
  if (F[S][T] != P_NOT_ASSIGNED) {
    return F[S][T];
  }
  constexpr auto MAX_N_PATHS = MAX_N_MEMBERS * MAX_N_MEMBERS;
  auto paths = StaticVector<edge_probability_t, MAX_N_PATHS>{};
  // Extracts all the direct paths in S -> Gv
  for (auto j : indices_in_bitset(N, S)) {
    for (auto i : gv.member_indices()) {
      paths.push_back(edge.p_cross[j][i]);
    }
  }
  auto p_success_in_this_turn = at_least_1_probability_r(paths);
  constexpr auto msg_pattern = "F[{1:#0{0}b}][{2:#0{0}b}]: Probability of this turn = {3:.4f}";
  // 2 : extra width of the prefix "0b"
  MYLOG_FMT_TRACE(msg_pattern, N + 2, T, S, p_success_in_this_turn);

  auto sum = 0.0_ep;
  auto to_p_internal = views::transform([&](auto j_pair) {
    auto [j0, j1] = j_pair;
    return gu.p_internal[j0][j1];
  });
  // Traverses all the non-empty subset of U - T
  for (auto S_next = U ^ T; S_next != 0; S_next = (S_next - 1) & (U ^ T)) {
    auto others = U ^ T ^ S_next;
    // For each dest in S_next, at least 1 path in S -> {dest} must pass
    auto p_S_to_S_next = 1.0_ep;
    for (auto dest : indices_in_bitset(N, S_next)) {
      p_S_to_S_next *= at_least_1_probability_r( //
          indices_in_bitset(N, S) | TRANSFORM_VIEW(gu.p_internal[_1][dest]));
    }
    // Each direct path in S -> others must fail
    auto p_S_to_others = at_least_1_probability_r(
        views::cartesian_product(indices_in_bitset(N, S), indices_in_bitset(N, others)) | to_p_internal);
    // Message propagates to S_next in the next turn
    auto p_next = get_merged_p_precise_impl(edge, gu, gv, F, T | S_next, S_next);
    constexpr auto msg_pattern = "T = {1:#0{0}b}, S = {2:#0{0}b}, S_next = {3:#0{0}b}, others = {4:#0{0}b}, "
                                 "p_S_to_S_next = {5:.4f}, p_S_to_others = {6:.4f}, p_next = {7:.4f}";
    // 2 : Extra width of the prefix "0b"
    MYLOG_FMT_TRACE(msg_pattern, N + 2, T, S, S_next, others, p_S_to_S_next, p_S_to_others, p_next);
    sum += p_S_to_S_next * (1.0_ep - p_S_to_others) * p_next;
  }
  BOOST_ASSERT_MSG(sum <= 1.0, "Implementation error: probability sum is out of range [0, 1].");
  // For memorized searching
  return F[T][S] = p_success_in_this_turn + (1.0_ep - p_success_in_this_turn) * sum;
}

auto get_merged_p_precise(const CoarsenedEdgeDetails& edge, const CoarsenedVertexDetails& gu,
                          const CoarsenedVertexDetails& gv, size_t j_bitset) -> edge_probability_t {
  auto F = make_initial_precise_state_container();
  return get_merged_p_precise_impl(edge, gu, gv, F, j_bitset, j_bitset);
}

auto get_merged_p_precise(const CoarsenedEdgeDetails& edge, const CoarsenedVertexDetails& gu,
                          const CoarsenedVertexDetails& gv, PreciseStateContainer& F, size_t j_bitset)
    -> edge_probability_t {
  return get_merged_p_precise_impl(edge, gu, gv, F, j_bitset, j_bitset);
}

auto merge_coarsened_wim_graph_edge(CoarsenedEdgeDetails& dest, const CoarsenedEdgeDetails* inv_dest,
                                    const CoarsenedVertexDetails& gu, const CoarsenedVertexDetails& gv,
                                    const CoarseningParams& params) -> void {
  auto N = gu.n_members();
  // Special case: Gu has only 1 member {u}
  if (N == 1) {
    auto Mv = gv.n_members();
    auto p_in_is_zero = get_normalized_heuristic_p_in(dest, inv_dest, gu, gv, params.in_out_heuristic_rule).is_all_zero;
    // If no in-edge to u, or all the in-edges of u are from Gv,
    // then there's no coarsened non-seed path accessible from Gu to Gv
    dest.merged.p = p_in_is_zero ? 0.0_ep : at_least_1_probability_r(dest.p_cross[0] | views::take(Mv));
    // However, the seed path is always enabled, whether or not there's in-edge to u.
    dest.merged.p_seed = at_least_1_probability_r(dest.p_seed_cross[0] | views::take(Mv));
    return;
  }
  BOOST_ASSERT_MSG(N == 2 || N == 3, "Group size other than 1, 2, 3 is not supported.");
  // Part (1): merged.p
  auto [p_in, p_in_all_zero] = get_normalized_heuristic_p_in(dest, inv_dest, gu, gv, params.in_out_heuristic_rule);
  MYLOG_FMT_TRACE("normalized in-heuristic of group {} = {::.4f}", gu.members, p_in);
  if (p_in_all_zero) {
    // Cannot serve as non-seed coarsened vertex since the group is isolated
    // (with no in-edge from other coarsened vertices)
    dest.merged.p = 0.0;
  } else if (params.edge_weight_rule != EdgeWeightRule::MERGED_PRECISE) {
    // SEPARATE_SIMPLE, MERGED_SIMPLE; SEPARATE_PRECISE: Calculates p_access(u) separately for each u in Gu
    // where p_access(u) = The (estimated) probability that message propagates from u to at least 1 vertex in Gv
    auto p_paths = StaticVector<edge_probability_t, MAX_N_MEMBERS>(N);
    if (params.edge_weight_rule == EdgeWeightRule::SEPARATE_PRECISE) {
      auto F = make_initial_precise_state_container();
      for (auto j0 : range(N)) {
        p_paths[j0] = get_merged_p_precise(dest, gu, gv, F, make_bitset_from_indices({j0}));
      }
    } else {
      for (auto j0 : range(N)) {
        p_paths[j0] = get_merged_p_simple(dest, gu, gv, j0, false);
      }
    }
    MYLOG_FMT_TRACE("p_paths from group {} to group {} = {::.4f}", gu.members, gv.members, p_paths);

    auto to_combined = TRANSFORM_VIEW(p_in[_1] * p_paths[_1]);
    if (params.edge_weight_rule == EdgeWeightRule::MERGED_SIMPLE) {
      // Conditional probability, where VS is a virtual source
      // with out-edge VS -> gu.members[i] with p = p_in[i] for each i = 0 ... gu.n_members() - 1
      dest.merged.p = at_least_1_probability_r(range(N) | to_combined) // P(VS -> gu.members -> gv.members)
                      / at_least_1_probability_r(p_in);                // P(VS -> gu.members)
    } else {
      BOOST_ASSERT(params.edge_weight_rule == EdgeWeightRule::SEPARATE_SIMPLE ||
                   params.edge_weight_rule == EdgeWeightRule::SEPARATE_PRECISE);
      // Weighted sum with p_in[i] as weight of gu.members[i]
      dest.merged.p = accumulate_sum(range(N) | to_combined);
    }
  } else {
    // MERGED_PRECISE: The probability that message propagates from VS to at least 1 vertex in Gv (via Gu)
    BOOST_ASSERT(params.edge_weight_rule == EdgeWeightRule::MERGED_PRECISE);
    auto U = (1zu << N) - 1; // Universal set {0 ... N-1}
    auto to_p_in = TRANSFORM_VIEW(p_in[_1]);
    auto to_p_not_in = TRANSFORM_VIEW(1.0_ep - p_in[_1]);
    // States are shared for the same pair (Gu, Gv)
    auto F = make_initial_precise_state_container();
    auto sum = 0.0_ep;
    // Traverses all the non-empty subsets of Gu
    for (auto S = 1zu; S <= U; S++) {
      // VS ---> u for each u in S
      auto p_enter_1 = accumulate_product(indices_in_bitset(N, S) | to_p_in);
      // VS -/-> u for each u in U - S
      auto p_enter_0 = accumulate_product(indices_in_bitset(N, U ^ S) | to_p_not_in);
      // S ---> Gv
      auto p_pass = get_merged_p_precise(dest, gu, gv, F, S);
      constexpr auto msg_pattern = "S = {1:#0{0}b}:\n\tp(VS ---> S) = {2:.4f}"
                                   "\n\tp(VS -/-> S) = {3:.4f}\n\tp(S ---> Gv) = {4:.4f}";
      MYLOG_FMT_TRACE(msg_pattern, N + 2, S, p_enter_1, p_enter_0, p_pass); // 2 : extra width of the prefix "0b"
      sum += p_enter_1 * p_enter_0 * p_pass;
    }
    BOOST_ASSERT_MSG(sum <= 1.0, "Implementation error: probability sum is out of range [0, 1].");
    // Conditional probability: P(VS -> gu.members -> gv.members) / P(VS -> gu.members)
    dest.merged.p = sum / at_least_1_probability_r(p_in);
  }

  // Part (2): merged.p_seed, which only considers the best candidate inside the group gu
  dest.merged.p_seed = [&]() {
    if (params.edge_seed_weight_rule == EdgeSeedWeightRule::BEST_SEED_INDEX) {
      auto res = get_merged_p_simple(dest, gu, gv, gu.best_seed_index, true);
      MYLOG_FMT_TRACE("p_seed from group {} to {} = {:.4f}, starting from the locally best seed candidate #{}", //
                      gu.members, gv.members, res, gu.members[gu.best_seed_index]);
      return res;
    }
    auto estimated_gain_as_seed = StaticVector<edge_probability_t, MAX_N_MEMBERS>(N);
    for (auto j0 : range(N)) {
      estimated_gain_as_seed[j0] = get_merged_p_simple(dest, gu, gv, j0, true);
    }
    MYLOG_FMT_TRACE("Estimated gain as seed of group {} = {}", gu.members, estimated_gain_as_seed);
    if (params.edge_seed_weight_rule == EdgeSeedWeightRule::AVERAGE) {
      return accumulate_sum(estimated_gain_as_seed) / N;
    }
    BOOST_ASSERT(params.edge_seed_weight_rule == EdgeSeedWeightRule::MAX);
    return ranges::max(estimated_gain_as_seed);
  }();

  // Part (3): Ensures 0.0 <= p <= p_seed <= 1.0
  dest.merged.p = std::clamp(dest.merged.p, 0.0_ep, 1.0_ep);
  dest.merged.p_seed = std::clamp(dest.merged.p_seed, dest.merged.p, 1.0_ep);
}

// If is_inverse, then the edge source -> target is mapped to Gv -> Gu
// Otherwise, the edge source -> target is mapped to Gu -> Gv
struct EdgeCoarseningItemInfo {
  bool is_inverse;
  vertex_id_t source;
  vertex_id_t target;
  WIMEdge w;

  auto dump() const {
    constexpr auto msg_pattern = "{{.is_inverse = {}, .source = {}, .target = {}, "
                                 ".w.p = {:.4f}, .w.p_seed = {:.4f}}}";
    return fmt::format(msg_pattern, is_inverse, source, target, w.p, w.p_seed);
  }
};

auto get_coarsened_wim_graph_edges(AdjacencyList<WIMEdge>& graph, const CoarseningDetails& details,
                                   const CoarseningParams& params) -> DirectedEdgeList<WIMEdge> {
  using ItemInfo = EdgeCoarseningItemInfo;
  using ItemType = VertexPairWithValue<const ItemInfo*>;

  auto m = graph.num_edges();
  auto item_info_values = make_reserved_vector<ItemInfo>(m);
  auto items = make_reserved_vector<ItemType>(m);

  for (auto [u, v, w] : graph::make_edge_range<0>(graph)) {
    auto [gu, gv] = details.to_group_id(u, v);
    if (gu == gv) {
      continue; // Skips the in-group edges
    }
    if (gu < gv) {
      auto pos = &item_info_values.emplace_back(ItemInfo{.is_inverse = false, .source = u, .target = v, .w = w});
      items.push_back(ItemType{.u = gu, .v = gv, .value = pos});
    } else {
      auto pos = &item_info_values.emplace_back(ItemInfo{.is_inverse = true, .source = u, .target = v, .w = w});
      items.push_back(ItemType{.u = gv, .v = gu, .value = pos});
    }
  }
  radix_sort_vertex_pairs(items, details.n_coarsened);
  MYLOG_TRACE([&]() {
    auto str = "Items after radix sorting: ["s;
    for (const auto& item : items) {
      str += fmt::format("\n\t{{.u = {}, .v = {}, .value -> {}}}", item.u, item.v, item.value->dump());
    }
    str += "\n]";
    return str;
  }());

  auto res = DirectedEdgeList<WIMEdge>{details.n_coarsened}; // Keeps coarsened |V| unchanged
  res.open_for_push_back();
  // Each chunk represents all the edges of Gu -> Gv (is_inverse == false) and Gv -> Gu (is_inverse == true)
  for (auto chunk : items | CHUNK_BY_VIEW(_1.u == _2.u && _1.v == _2.v)) {
    auto [gu, gv] = chunk[0].vertices();
    auto [Mu, Mv] = details.group_id_to_size(gu, gv);
    auto E_uv = CoarsenedEdgeDetails{.n_members_left = Mu, .n_members_right = Mv}; // Gu -> Gv
    auto E_vu = CoarsenedEdgeDetails{.n_members_left = Mv, .n_members_right = Mu}; // Gv -> Gu

    auto has_uv = false;
    auto has_vu = false;
    for (auto info : chunk | views::transform(&ItemType::value)) {
      auto [j_source, j_target] = details.to_index_in_group(info->source, info->target);
      if (info->is_inverse) {
        // The edge source -> target is mapped to Gv -> Gu
        E_vu.p_cross[j_source][j_target] = info->w.p;
        E_vu.p_seed_cross[j_source][j_target] = info->w.p_seed;
        has_vu = true;
      } else {
        // The edge source -> target is mapped to Gu -> Gv
        E_uv.p_cross[j_source][j_target] = info->w.p;
        E_uv.p_seed_cross[j_source][j_target] = info->w.p_seed;
        has_uv = true;
      }
    }
    if (has_uv) {
      merge_coarsened_wim_graph_edge(E_uv, &E_vu, details.groups[gu], details.groups[gv], params);
      res.push_back(gu, gv, E_uv.merged);
      MYLOG_FMT_TRACE("Details of coarsened edge {} -> {}: {:4}", gu, gv, E_uv);
    }
    if (has_vu) {
      merge_coarsened_wim_graph_edge(E_vu, &E_uv, details.groups[gv], details.groups[gu], params);
      res.push_back(gv, gu, E_vu.merged);
      MYLOG_FMT_TRACE("Details of coarsened edge {} -> {}: {:4}", gv, gu, E_vu);
    }
  }
  res.close_for_push_back();
  return res;
}

auto get_coarsened_wim_graph_vertex_weights(const CoarseningDetails& details,
                                            std::span<const vertex_weight_t> vertex_weights,
                                            const CoarseningParams& params) -> std::vector<vertex_weight_t> {
  auto coarsened_vertex_weights = std::vector<vertex_weight_t>(details.n_coarsened, 0.0);
  // vc = Coarsened vertex index in range [0, Nc)
  for (auto [vc, gc] : details.groups | views::enumerate) {
    auto [p_in, p_in_all_zero] = get_normalized_heuristic_p_in(gc, params.in_out_heuristic_rule);
    auto& dest = coarsened_vertex_weights[vc];

    if (params.vertex_weight_rule == VertexWeightRule::SUM) {
      // SUM: Simply sums all the vertex weights up
      dest = accumulate_sum(gc.members | TRANSFORM_VIEW(vertex_weights[_1]));
    } else if (params.vertex_weight_rule == VertexWeightRule::AVERAGE || p_in_all_zero) {
      // (1) AVERAGE: Takes arithmetic average
      // (2) Special case for AVERAGE and AVERAGE_BY_PATHS: Takes arithmetic average if no in-heuristic
      dest = accumulate_sum(gc.members | TRANSFORM_VIEW(vertex_weights[_1])) / gc.n_members();
    } else if (gc.members.size() == 1) {
      // AVERAGE_BY_PATHS for group of size 1
      dest = vertex_weights[gc.members[0]];
    } else if (gc.members.size() == 2) {
      // AVERAGE_BY_PATHS for group of size 2
      dest = 0.0_vw;
      for (auto [i1, v1] : gc.members | views::enumerate) {
        // (1) v1; (2) v0 -> v1
        auto i0 = get_next_j<2>(i1);
        auto paths = {p_in[i1], p_in[i0] * gc.p_internal[i0][i1]};
        dest += vertex_weights[v1] * at_least_1_probability_r(paths);
      }
    } else {
      // AVERAGE_BY_PATHS for group of size 3
      dest = 0.0_vw;
      BOOST_ASSERT(gc.members.size() == 3);
      for (auto [i2, v2] : gc.members | views::enumerate) {
        auto [i0, i1] = get_next_j<3>(i2);
        auto prob = at_least_1_probability(                            //
            p_in[i2],                                                  // (1) v2
            p_in[i0] * gc.p_internal[i0][i2],                          // (2) v0 -> v2
            p_in[i0] * gc.p_internal[i0][i1] * gc.p_internal[i1][i2],  // (3) v0 -> v1 -> v2
            p_in[i1] * gc.p_internal[i1][i2],                          // (4) v1 -> v2
            p_in[i1] * gc.p_internal[i1][i0] * gc.p_internal[i0][i2]); // (5) v1 -> v0 -> v2
        dest += vertex_weights[v2] * prob;
      }
    }
  }
  return coarsened_vertex_weights;
}

template <class ResultType>
auto coarsen_wim_graph_w_impl(AdjacencyList<WIMEdge>& graph, InvAdjacencyList<WIMEdge>& inv_graph,
                              std::span<const vertex_weight_t> vertex_weights,             //
                              vertex_id_t n_groups, std::span<const vertex_id_t> group_id, //
                              const CoarseningParams& params) noexcept -> ResultType {
  MYLOG_FMT_DEBUG("Parameters: {:4}", params);
  auto n = graph::num_vertices(graph);
  // d_in: unit weight by default
  auto details = CoarseningDetails{.n = n,
                                   .n_coarsened = n_groups,
                                   .group_id = std::vector(group_id.begin(), group_id.end()),
                                   .index_in_group = std::vector<vertex_id_t>(n),
                                   .groups = std::vector<CoarsenedVertexDetails>(n_groups)};
  MYLOG_FMT_DEBUG("Starts coarsening: n = {}, n_groups = {}", details.n, details.n_coarsened);

  auto total_time_used = 0.0;
  auto timer = nw::util::seconds_timer{};
  auto step_block = [&](int step_index, std::string_view step_description, auto&& step_fn) {
    timer.start();
    std::invoke(step_fn);
    timer.stop();
    total_time_used += timer.elapsed();
    constexpr auto msg_pattern = "Done coarsening step {} after Mongoose matching: {}."
                                 "\n\tTime usage: {:.3f} in current step, {:.3f} in total.";
    ELOGFMT(INFO, msg_pattern, step_index, step_description, timer.elapsed(), total_time_used);
  };

  // Step 1: groups[].members & groups[].vertex_weights & groups[].heuristics_in & index_in_group
  // Vertex v is placed in group g
  step_block(1, "assigning group members", [&]() {
    for (auto [v, g] : views::enumerate(group_id)) {
      auto& cur_group = details.groups[g];
      details.index_in_group[v] = static_cast<vertex_id_t>(cur_group.members.size());
      cur_group.members.push_back(v);
      cur_group.vertex_weights.push_back(vertex_weights[v]);
      // Heuristics
      cur_group.heuristics_in.push_back( //
          get_heuristic_in(inv_graph, group_id, vertex_weights, v, params.in_out_heuristic_rule));
      cur_group.heuristics_out.push_back( //
          get_heuristic_out(graph, group_id, vertex_weights, v, params.in_out_heuristic_rule));
      cur_group.heuristics_out_seed.push_back( //
          get_heuristic_out_seed(graph, group_id, vertex_weights, v, params.in_out_heuristic_rule));
    }
  });
  // Step 2: groups[].p_internal
  step_block(2, "assigning p_internal of groups", [&]() {
    for (auto [u, v, w] : graph::make_edge_range<0>(graph)) {
      auto g = group_id[u];
      if (g == group_id[v]) {
        auto [ju, jv] = details.to_index_in_group(u, v);
        details.groups[g].p_internal[ju][jv] = w.p;
        details.groups[g].p_seed_internal[ju][jv] = w.p_seed;
      }
    }
  });
  // Step 3: Calculates the best seed locally in each group
  step_block(3, "selecting best seed candidate in each group", [&]() {
    for (auto& g : details.groups) {
      g.best_seed_index = select_best_seed_in_group(g, params.in_out_heuristic_rule).index_in_group;
    }
  });
  // Step 4: Generates coarsened vertex weights
  auto coarsened_vertex_weights = std::vector<vertex_weight_t>{};
  step_block(4, "coarsening vertex weights", [&]() {
    coarsened_vertex_weights = get_coarsened_wim_graph_vertex_weights(details, vertex_weights, params);
  });
  // Step 5: Finally, builds the result from all the components
  auto coarsened_edge_list = DirectedEdgeList<WIMEdge>{};
  step_block(5, "building coarsened edge list", [&]() {
    coarsened_edge_list = get_coarsened_wim_graph_edges(graph, details, params);
    BOOST_ASSERT_MSG(graph::num_vertices(coarsened_edge_list) == details.n_coarsened,
                     "Mismatch between expected and actual # of coarsened vertices.");
  });

  auto res_details = [&]() {
    if constexpr (std::is_same_v<ResultType, CoarsenGraphDetailedResult>) {
      // Returns the details directly
      return std::move(details);
    } else if constexpr (std::is_same_v<ResultType, CoarsenGraphBriefResult>) {
      // Transforms detailed to brief information
      auto brief = CoarseningBrief{
          .n = n, .n_coarsened = n_groups, .groups = make_reserved_vector<CoarsenedVertexBrief>(n_groups)};
      for (const auto& u : details.groups) {
        brief.groups.push_back({.members = u.members, .best_seed_index = u.best_seed_index});
      }
      return brief;
    } else {
      static_assert(rfl::always_false_v<ResultType>, "Invalid result type.");
    }
  }();
  auto coarsened = WIMAdjacencyListPair{
      .adj_list = AdjacencyList<WIMEdge>{coarsened_edge_list},
      .inv_adj_list = InvAdjacencyList<WIMEdge>{coarsened_edge_list},
      .vertex_weights = std::move(coarsened_vertex_weights),
  };
  return ResultType{
      .coarsened = std::move(coarsened),
      .details = std::move(res_details),
  };
}

template <class ResultType>
auto coarsen_wim_graph_w_framework(const AdjacencyList<WIMEdge>& graph, const InvAdjacencyList<WIMEdge>& inv_graph,
                                   std::span<const vertex_weight_t> vertex_weights,
                                   const CoarseningParams& params) noexcept -> ResultType {
  // Step 1:
  MYLOG_FMT_DEBUG("Coarsening: |V|, |E| of graph = {}; |V|, |E| of inv_graph = {}; Size of vertex_weights = {}", //
                  graph_n_m(graph), graph_n_m(inv_graph), vertex_weights.size());
  auto bidir_graph = merge_wim_edge_to_undirected(graph, params);
  MYLOG_FMT_DEBUG("Coarsening: |V|, |E| of bidir_graph = {}", graph_n_m(bidir_graph));
  // Step 2
  auto [n_groups, group_id] = mongoose_match(bidir_graph, params);
  MYLOG_FMT_DEBUG("Coarsening: n_groups = {}; Size of group_id = {}", n_groups, group_id.size());
  // Step 3
  return coarsen_wim_graph_w_impl<ResultType>(as_non_const(graph), as_non_const(inv_graph), //
                                              vertex_weights, n_groups, group_id, params);
}

template <class CoarsenResultType>
auto expand_wim_seed_vertices_w_impl(const AdjacencyList<WIMEdge>& graph,
                                     std::span<const vertex_weight_t> vertex_weights,
                                     const CoarsenResultType& coarsening_result,
                                     std::span<const vertex_id_t> coarsened_seeds,
                                     const ExpandingParams& params) noexcept -> rfl::Result<ExpandSeedResult> {
  auto n = coarsening_result.details.n;            // # of vertices BEFORE coarsening
  auto nc = coarsening_result.details.n_coarsened; // # of vertices AFTER coarsening
  BOOST_ASSERT_MSG(ranges::size(coarsened_seeds) <= nc, "Expects more than nc vertices as seeds");
  BOOST_ASSERT_MSG(ranges::max(coarsened_seeds) < nc, "Indices of seeds are out of range [0, nc).");

  auto expanded_seeds = make_reserved_vector<vertex_id_t>(coarsened_seeds.size());
  auto make_expanded_seeds_as_s_local = [&]() {
    for (auto sc : coarsened_seeds) {
      auto& group = coarsening_result.details.groups[sc];
      expanded_seeds.push_back(group.members[group.best_seed_index]);
    }
  };
  MYLOG_FMT_DEBUG("coarsened_seeds = {}", coarsened_seeds);

  if (params.seed_expanding_rule == SeedExpandingRule::LOCAL) {
    // S_LOCAL: Only considers information inside current group
    make_expanded_seeds_as_s_local();
  } else if (params.seed_expanding_rule == SeedExpandingRule::SIMULATIVE) {
    // S_SIMULATIVE: Picks one-by-one separately by simulation result of single expanded seed vertex
    for (auto sc : coarsened_seeds) {
      auto& group = coarsening_result.details.groups[sc];
      auto estimated_gain = StaticVector<double, MAX_N_MEMBERS>{};
      for (auto v : group.members) {
        auto sim_res = wim_simulate_w(graph, vertex_weights, {n, {v}}, **params.simulation_try_count);
        RFL_RETURN_ON_ERROR(sim_res);
        estimated_gain.push_back(*sim_res);
      }
      auto best_expanded_seed = group.members[ranges::max_element(estimated_gain) - estimated_gain.begin()];
      constexpr auto msg_pattern = "Expanding seeds: Estimated gain of group {} "
                                   "(coarsened to vertex #{}) = {::.4f} (best member: {})";
      MYLOG_FMT_DEBUG(msg_pattern, group.members, sc, estimated_gain, best_expanded_seed);
      expanded_seeds.push_back(best_expanded_seed);
    }
  } else {
    // S_ITERATIVE: Result of S_SEPARATE as initial solution
    BOOST_ASSERT(params.seed_expanding_rule == SeedExpandingRule::ITERATIVE);
    make_expanded_seeds_as_s_local();
    auto sim_res = simulate(graph, vertex_weights, expanded_seeds, nullptr, **params.simulation_try_count);
    RFL_RETURN_ON_ERROR(sim_res);
    MYLOG_FMT_DEBUG("Initial expanded seed set: {}\n\tWith simulation result = {:.4f}", expanded_seeds, *sim_res);
    // Changable only of at least 2 vertices in the group that is expanded from the coarsened seed vertex.
    auto changable_group_indices = std::vector<size_t>{};
    for (auto [i, sc] : coarsened_seeds | views::enumerate) {
      if (coarsening_result.details.groups[sc].members.size() >= 2) {
        changable_group_indices.push_back(i);
      }
    }
    // Special case: each coarsened seed is expanded to only 1 vertex.
    if (changable_group_indices.empty()) {
      goto expand_wim_seed_vertices_w_return;
    }
    // Each iteration randomly changes one seed vertex
    for (auto iteration_index : range(**params.n_iterations)) {
      // Ensures that each group is attempted if # of iterations >= # of groups corresponding to coarsened seeds
      if (iteration_index % changable_group_indices.size() == 0) {
        ranges::shuffle(changable_group_indices, rand_engine);
      }
      auto si = changable_group_indices[iteration_index % changable_group_indices.size()];
      auto& group = coarsening_result.details.groups[coarsened_seeds[si]];
      // Ensures next != current expanded seed
      auto [prev, next] = [&]() {
        auto res = rand_element(group.members);
        while (res == expanded_seeds[si]) {
          res = rand_element(group.members);
        }
        return std::tuple{expanded_seeds[si], res};
      }();
      MYLOG_FMT_DEBUG("Iteration #{}: Attempting to change {} to {} (in group {} coarsened to vertex #{})",
                      iteration_index, prev, next, group.members, coarsened_seeds[si]);
      expanded_seeds[si] = next;
      auto next_sim_res = simulate(graph, vertex_weights, expanded_seeds, nullptr, **params.simulation_try_count);
      RFL_RETURN_ON_ERROR(next_sim_res);
      MYLOG_FMT_DEBUG("Iteration #{}: Candidate {} vs. {}: Simulation result: {:.4f} vs. {:.4f}", iteration_index, next,
                      prev, *next_sim_res, *sim_res);
      if (*next_sim_res > *sim_res) {
        sim_res = next_sim_res; // Changes prev to next
      } else {
        expanded_seeds[si] = prev; // Otherwise, restores current expanded seed
      }
      MYLOG_FMT_DEBUG("Expanded seed set after iteration #{}: {}", iteration_index, expanded_seeds);
    }
  }
expand_wim_seed_vertices_w_return:
  return ExpandSeedResult{.expanded_seeds = std::move(expanded_seeds)};
}
} // namespace

// ---- Step 1 ----

auto merge_wim_edge_to_undirected(const AdjacencyList<WIMEdge>& graph, const CoarseningParams& params) noexcept
    -> AdjacencyList<edge_probability_t> {
  // Workaround: const reference triggers a bug in nwgraph. The graph is not changed actually.
  return merge_edge_to_undirected_generic(const_cast<AdjacencyList<WIMEdge>&>(graph), params);
}

auto merge_wbim_edge_to_undirected(const AdjacencyList<WBIMEdge>& graph, const CoarseningParams& params) noexcept
    -> AdjacencyList<edge_probability_t> {
  // Workaround: const reference triggers a bug in nwgraph. The graph is not changed actually.
  return merge_edge_to_undirected_generic(const_cast<AdjacencyList<WBIMEdge>&>(graph), params);
}

// ---- Step 2 ----

auto mongoose_match(const AdjacencyList<edge_probability_t>& graph, const CoarseningParams& params) noexcept
    -> MongooseMatchResult {
  return mongoose_match_with_seeds(graph, {}, params);
}

auto mongoose_match_with_seeds(const AdjacencyList<edge_probability_t>& graph, std::span<const vertex_id_t> seeds,
                               const CoarseningParams& params) noexcept -> MongooseMatchResult {
  auto n = graph::num_vertices(graph);
  auto seq = std::vector<vertex_id_t>(n);
  std::ranges::iota(seq, 0);
  MYLOG_FMT_DEBUG("Starts Mongoose matching algorithm. |V| = {}", n);

  auto match = std::vector<vertex_id_t>(n, -1);
  auto adopt = std::vector<vertex_id_t>(n, -1);
  auto make_match = [&](vertex_id_t u, vertex_id_t v) {
    match[u] = v;
    match[v] = u;
  };

  // Step 1: Neighbor matching
  ranges::shuffle(seq, rand_engine);
  for (auto s : seeds) {
    BOOST_ASSERT_MSG(s >= 0 && s < n, "Seed vertex out of range [0, n).");
    match[s] = s; // Marks the seeds as single-matched
  }
  for (auto v : seq) {
    if (match[v] != -1) {
      continue; // Skips the matched vertices (including seeds)
    }
    auto best_p = p_infinity_by_rule(params.neighbor_match_rule);
    auto best_u = -1_vid;
    for (auto [u, p] : graph[v]) {
      // Selects an unmatched neighbor with the best edge weight
      if (match[u] == -1 && is_better_p(best_p, p, params.neighbor_match_rule)) {
        best_p = p;
        best_u = u;
      }
    }
    if (best_u == -1) {
      MYLOG_FMT_TRACE("#{} is isolated during neighbor matching.", v);
      match[v] = v;
    } else {
      MYLOG_FMT_TRACE("#{} matches with #{}", v, best_u);
      make_match(v, best_u);
    }
  }
  MYLOG_FMT_TRACE("Step 1 done: Neighbor matching.\n{}", DUMP_INDEX_ARRAY(match));

  // Step 2: Brotherly matching & Adoptive matching
  ranges::shuffle(seq, rand_engine);
  for (auto s : seeds) {
    BOOST_ASSERT_MSG(s >= 0 && s < n, "Seed vertex out of range [0, n).");
    match[s] = -1; // Uses -1 to mark seeds
  }
  for (auto v : seq) {
    if (match[v] != v) {
      continue; // Skips the matched vertices or seeds
    }
    auto v_pivot = [&] {
      auto best_p = p_infinity_by_rule(params.neighbor_match_rule);
      auto best_u = -1_vid;
      for (auto [u, p] : graph[v]) {
        BOOST_ASSERT_MSG(match[u] == -1 || match[u] != u, "Unexpected non-maximal matching in step 1.");
        if (match[u] != -1 && is_better_p(best_p, p, params.neighbor_match_rule)) {
          best_u = u;
          best_p = p;
        }
      }
      return best_u;
    }();
    MYLOG_FMT_TRACE("v = {}: pivot = {}", v, v_pivot);
    if (v_pivot == -1) {
      MYLOG_FMT_TRACE("Found isolated vertex {} during brotherly matching.", v);
      continue; // Special case when v is an isloated vertex, or all the neighbors are seeds
    }
    // Step 2.1 for each vertex v: Neighbor matching
    auto u_last = -1_vid;
    for (auto u : graph[v_pivot] | views::keys | views::filter(LAMBDA_1(match[_1] == _1))) {
      if (u_last != -1) {
        MYLOG_FMT_TRACE("Brotherly matching around pivot #{}: {} -- {}", v_pivot, u_last, u);
        make_match(u_last, u);
        u_last = -1;
      } else {
        u_last = u;
      }
    }
    if (u_last == -1) {
      continue;
    }
    // Step 2.2 for each vertex v: Adoptive matching
    auto x = match[v_pivot];
    MYLOG_FMT_TRACE("Begins adoptive matching for #{}, with pivot = {}, match[pivot] = {}", u_last, v_pivot, x);
    // v1 -> v2 : v2 adopts v1;  v1 -- v2 : (v1, v2) forms a match
    if (adopt[v_pivot] != -1) {
      constexpr auto msg_pattern_1 = "Adoptive matching (1): {0} -> {1} -- {2}  ==>  {0} -- {1}, {2} -- {3}";
      MYLOG_FMT_TRACE(msg_pattern_1, adopt[v_pivot], v_pivot, x, u_last);
      make_match(u_last, x);
      make_match(adopt[v_pivot], v_pivot);
      adopt[v_pivot] = -1;
    } else if (adopt[x] != -1) {
      constexpr auto msg_pattern_2 = "Adoptive matching (2): {0} -> {1} -- {2}  ==>  {0} -- {1}, {2} -- {3}";
      MYLOG_FMT_TRACE(msg_pattern_2, adopt[x], x, v_pivot, u_last);
      make_match(u_last, v_pivot);
      make_match(adopt[x], x);
      adopt[x] = -1;
    } else {
      constexpr auto msg_pattern_3 = "Adoptive matching (3): {0} -- {1}  ==>  {0} -- {1} <- {2}";
      MYLOG_FMT_TRACE(msg_pattern_3, x, v_pivot, u_last);
      match[u_last] = v_pivot;
      adopt[v_pivot] = u_last;
    }
  }
  MYLOG_TRACE("Step 2 done: Brotherly matching & Adoptive matching.\n{}\n{}", //
              DUMP_INDEX_ARRAY(match), DUMP_INDEX_ARRAY(adopt));

  // Step 3: Groups by matching result
  auto res = MongooseMatchResult{.n_groups = 0, .group_id = std::vector<vertex_id_t>(n, -1)};
  for (auto v : range(n)) {
    if (res.group_id[v] != -1) {
      continue; // Ignores the already grouped vertices
    }
    if (match[v] == -1 || match[v] == v) {
      // (1) v is isolated, or v is a seed vertex
      MYLOG_FMT_TRACE("Isolated vertex found during grouping: {}", v);
      res.group_id[v] = res.n_groups;
    } else {
      auto x = match[v];
      if (match[x] == v) {
        // (2) ? -> v -- x <- ?
        res.group_id[v] = res.group_id[x] = res.n_groups;
        if (adopt[v] != -1) {
          res.group_id[adopt[v]] = res.n_groups;
        } else if (adopt[x] != -1) {
          res.group_id[adopt[x]] = res.n_groups;
        }
      } else {
        // (3) v -> x -- ?
        BOOST_ASSERT_MSG(match[x] != x && match[match[x]] == x,
                         "Unexpected behavior: (x, match[x]) is not matched bidirectionally.");
        res.group_id[v] = res.group_id[x] = res.group_id[match[x]] = res.n_groups;
      }
    }
    res.n_groups += 1; // Creates a new group
  }
  return res;
}

// ---- Step 3 ----

auto coarsen_wim_graph_with_match_result(const AdjacencyList<WIMEdge>& graph,
                                         const InvAdjacencyList<WIMEdge>& inv_graph, vertex_id_t n_groups,
                                         std::span<const vertex_id_t> group_id, //
                                         const CoarseningParams& params) noexcept -> CoarsenGraphBriefResult {
  auto vertex_weights = std::vector<vertex_weight_t>(graph::num_vertices(graph), 1.0);
  return coarsen_wim_graph_with_match_result_w(graph, inv_graph, vertex_weights, n_groups, group_id, params);
}

auto coarsen_wim_graph_with_match_result_w(const AdjacencyList<WIMEdge>& graph,
                                           const InvAdjacencyList<WIMEdge>& inv_graph,
                                           std::span<const vertex_weight_t> vertex_weights, vertex_id_t n_groups,
                                           std::span<const vertex_id_t> group_id,
                                           const CoarseningParams& params) noexcept -> CoarsenGraphBriefResult {
  return coarsen_wim_graph_w_impl<CoarsenGraphBriefResult>(const_cast<AdjacencyList<WIMEdge>&>(graph),
                                                           const_cast<InvAdjacencyList<WIMEdge>&>(inv_graph),
                                                           vertex_weights, n_groups, group_id, params);
}

auto coarsen_wim_graph_with_match_result_d(const AdjacencyList<WIMEdge>& graph,
                                           const InvAdjacencyList<WIMEdge>& inv_graph, vertex_id_t n_groups,
                                           std::span<const vertex_id_t> group_id,
                                           const CoarseningParams& params) noexcept -> CoarsenGraphDetailedResult {
  auto vertex_weights = std::vector<vertex_weight_t>(graph::num_vertices(graph), 1.0);
  return coarsen_wim_graph_with_match_result_d_w(graph, inv_graph, vertex_weights, n_groups, group_id, params);
}

auto coarsen_wim_graph_with_match_result_d_w(const AdjacencyList<WIMEdge>& graph,
                                             const InvAdjacencyList<WIMEdge>& inv_graph,
                                             std::span<const vertex_weight_t> vertex_weights, vertex_id_t n_groups,
                                             std::span<const vertex_id_t> group_id,
                                             const CoarseningParams& params) noexcept -> CoarsenGraphDetailedResult {
  return coarsen_wim_graph_w_impl<CoarsenGraphDetailedResult>(const_cast<AdjacencyList<WIMEdge>&>(graph),
                                                              const_cast<InvAdjacencyList<WIMEdge>&>(inv_graph),
                                                              vertex_weights, n_groups, group_id, params);
}

auto coarsen_wim_graph(const AdjacencyList<WIMEdge>& graph, const InvAdjacencyList<WIMEdge>& inv_graph,
                       const CoarseningParams& params) noexcept -> CoarsenGraphBriefResult {
  auto vertex_weights = std::vector<vertex_weight_t>(graph::num_vertices(graph), 1.0);
  return coarsen_wim_graph_w(graph, inv_graph, vertex_weights, params);
}

auto coarsen_wim_graph_w(const AdjacencyList<WIMEdge>& graph, const InvAdjacencyList<WIMEdge>& inv_graph,
                         std::span<const vertex_weight_t> vertex_weights, const CoarseningParams& params) noexcept
    -> CoarsenGraphBriefResult {
  return coarsen_wim_graph_w_framework<CoarsenGraphBriefResult>(graph, inv_graph, vertex_weights, params);
}

auto coarsen_wim_graph_d(const AdjacencyList<WIMEdge>& graph, const InvAdjacencyList<WIMEdge>& inv_graph,
                         const CoarseningParams& params) noexcept -> CoarsenGraphDetailedResult {
  auto vertex_weights = std::vector<vertex_weight_t>(graph::num_vertices(graph), 1.0);
  return coarsen_wim_graph_d_w(graph, inv_graph, vertex_weights, params);
}

auto coarsen_wim_graph_d_w(const AdjacencyList<WIMEdge>& graph, const InvAdjacencyList<WIMEdge>& inv_graph,
                           std::span<const vertex_weight_t> vertex_weights, const CoarseningParams& params) noexcept
    -> CoarsenGraphDetailedResult {
  return coarsen_wim_graph_w_framework<CoarsenGraphDetailedResult>(graph, inv_graph, vertex_weights, params);
}

// ---- Step 4 ----

auto select_best_seed_in_group(const CoarsenedVertexDetails& v, InOutHeuristicRule rule) noexcept
    -> SelectBestSeedResult {
  auto m = v.members.size();
  if (m == 1) {
    return {.index_in_group = 0};
  }
  BOOST_ASSERT(m <= 3);
  auto [p_out, p_out_all_zero] = get_normalized_heuristic_p_out(v, rule);
  auto [p_out_seed, p_out_seed_all_zero] = get_normalized_heuristic_p_out_seed(v, rule);
  auto to_estimated_gain = [&](size_t j0) {
    auto res = edge_probability_t{};
    if (v.n_members() == 2) {
      auto j1 = get_next_j<2>(j0);
      res = p_out_seed[j0] + v.p_seed_internal[j0][j1] * p_out[j1];
    } else {
      BOOST_ASSERT(v.n_members() == 3);
      auto [j1, j2] = get_next_j<3>(j0);
      auto p1 = at_least_1_probability(                      //
          v.p_seed_internal[j0][j1],                         // j0 -> j1
          v.p_seed_internal[j0][j2] * v.p_internal[j2][j1]); // j0 -> j2 -> j1
      auto p2 = at_least_1_probability(                      //
          v.p_seed_internal[j0][j2],                         // j0 -> j2
          v.p_seed_internal[j0][j1] * v.p_internal[j1][j2]); // j0 -> j1 -> j2
      res = p_out_seed[j0] + p1 * p_out[j1] + p2 * p_out[j2];
    }
    return res;
  };
  auto estimated_gain = [&]() {
    auto view = range(m) | views::transform(to_estimated_gain);
    return std::vector(view.begin(), view.end());
  }();
  MYLOG_FMT_TRACE("Best seed selection: estimated gain for group {} = {::.4f}", //
                  v.members, estimated_gain);
  return {.index_in_group = static_cast<size_t>(ranges::max_element(estimated_gain) - estimated_gain.begin())};
}

auto expand_wim_seed_vertices(const AdjacencyList<WIMEdge>& graph, const CoarsenGraphBriefResult& coarsening_result,
                              std::span<const vertex_id_t> coarsened_seeds, const ExpandingParams& params) noexcept
    -> rfl::Result<ExpandSeedResult> {
  auto vertex_weights = std::vector<vertex_weight_t>(graph::num_vertices(graph), 1.0);
  return expand_wim_seed_vertices_w(graph, vertex_weights, coarsening_result, coarsened_seeds, params);
}

auto expand_wim_seed_vertices_w(const AdjacencyList<WIMEdge>& graph, std::span<const vertex_weight_t> vertex_weights,
                                const CoarsenGraphBriefResult& coarsening_result,
                                std::span<const vertex_id_t> coarsened_seeds, const ExpandingParams& params) noexcept
    -> rfl::Result<ExpandSeedResult> {
  return expand_wim_seed_vertices_w_impl(graph, vertex_weights, coarsening_result, coarsened_seeds, params);
}

auto further_expand_wim_seed_vertices(const CoarsenGraphBriefResult& last_result,
                                      const CoarsenGraphBriefResult& cur_result,
                                      std::span<const vertex_id_t> coarsened_seeds,
                                      const ExpandingParams& params) noexcept -> rfl::Result<ExpandSeedResult> {
  return expand_wim_seed_vertices_w(last_result.coarsened.adj_list, last_result.coarsened.vertex_weights, cur_result,
                                    coarsened_seeds, params);
}

auto expand_wim_seed_vertices_d(const AdjacencyList<WIMEdge>& graph,
                                const CoarsenGraphDetailedResult& coarsening_result,
                                std::span<const vertex_id_t> coarsened_seeds, const ExpandingParams& params) noexcept
    -> rfl::Result<ExpandSeedResult> {
  auto vertex_weights = std::vector<vertex_weight_t>(graph::num_vertices(graph), 1.0);
  return expand_wim_seed_vertices_d_w(graph, vertex_weights, coarsening_result, coarsened_seeds, params);
}

auto expand_wim_seed_vertices_d_w(const AdjacencyList<WIMEdge>& graph, std::span<const vertex_weight_t> vertex_weights,
                                  const CoarsenGraphDetailedResult& coarsening_result,
                                  std::span<const vertex_id_t> coarsened_seeds, const ExpandingParams& params) noexcept
    -> rfl::Result<ExpandSeedResult> {
  return expand_wim_seed_vertices_w_impl(graph, vertex_weights, coarsening_result, coarsened_seeds, params);
}

auto further_expand_wim_seed_vertices_d(const CoarsenGraphDetailedResult& last_result,
                                        const CoarsenGraphDetailedResult& cur_result,
                                        std::span<const vertex_id_t> coarsened_seeds,
                                        const ExpandingParams& params) noexcept -> rfl::Result<ExpandSeedResult> {
  return expand_wim_seed_vertices_d_w(last_result.coarsened.adj_list, last_result.coarsened.vertex_weights, cur_result,
                                      coarsened_seeds, params);
}
