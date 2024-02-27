#include "coarsening.h"
#include "dump.h"
#include "utils/boost_assert.h"
#include "utils/graph.h"
#include "utils/static_vector.h"
#include "wim.h"
#include <fmt/ranges.h>
#include <map>
#include <numeric>
#include <nwgraph/adaptors/edge_range.hpp>
#include <ylt/easylog.hpp>

namespace {
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

auto initial_p(NeighborMatchRule rule) -> edge_probability_t {
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

constexpr auto P_NOT_ASSIGNED = -1.0_ep;

auto make_initial_coarsened_edge_details(size_t n_members_left, size_t n_members_right) -> CoarsenedEdgeDetails {
  return {.n_members_left = n_members_left,
          .n_members_right = n_members_right,
          .p_cross = {},
          .p_seed_cross = {},
          .merged = {.p = P_NOT_ASSIGNED, .p_seed = P_NOT_ASSIGNED}};
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

template <is_edge_property E>
auto merge_edge_to_undirected_generic(AdjacencyList<E>& graph, const CoarseningParams& params)
    -> AdjacencyList<edge_probability_t> {
  // Ensures |V| is unchanged
  auto edge_list = UndirectedEdgeList<edge_probability_t>{graph::num_vertices(graph)};
  edge_list.open_for_push_back();

  using VertexPair = std::pair<vertex_id_t, vertex_id_t>;
  auto edges_map = std::map<VertexPair, edge_probability_t>{};
  // Makes sure that each key (u, v) satisfied u < v
  for (auto [u, v, w] : graph::make_edge_range<0>(graph)) {
    auto it = (u < v) ? edges_map.find({u, v}) : edges_map.find({v, u});
    if (it == edges_map.end()) {
      edges_map.insert({(u < v) ? std::pair{u, v} : std::pair{v, u}, w.p});
    } else {
      it->second = merge_p(it->second, w.p, params.neighbor_match_rule);
    }
  }
  for (auto [v_pair, p] : edges_map) {
    edge_list.push_back(get<0>(v_pair), get<1>(v_pair), p);
  }
  edge_list.close_for_push_back();

  auto res = AdjacencyList<edge_probability_t>{edge_list};
  BOOST_ASSERT_MSG(graph::num_vertices(res) == graph::num_vertices(graph), "# of vertices mismatch.");
  return res;
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

struct GetHeuristicPInResult {
  CoarsenedVertexDetails::HeuristicsContainer p_in;
  bool is_all_zero;
};

auto get_heuristic_p_in(const CoarsenedEdgeDetails& edge, const CoarsenedEdgeDetails* inv_edge,
                        const CoarsenedVertexDetails& gu, const CoarsenedVertexDetails& gv, InOutHeuristicRule rule)
    -> GetHeuristicPInResult {
  auto res = GetHeuristicPInResult{.p_in = gu.heuristics_in, .is_all_zero = false};
  if (inv_edge != nullptr) {
    // Removes the directed edge gu.members[j] <-- gv.members[i]
    for (auto [j, i] : views::cartesian_product(gu.member_indices(), gv.member_indices())) {
      if (rule == InOutHeuristicRule::P) {
        res.p_in[j] -= inv_edge->p_cross[i][j];
      } else if (rule == InOutHeuristicRule::W) {
        res.p_in[j] -= inv_edge->p_cross[i][j] * gv.vertex_weights[i];
      } else {
        res.p_in[j] -= 1.0_ep;
      }
    }
    ELOGFMT(TRACE, "In-heuristics of group {} after removing back edges from {} = {}", //
            gu.members, gv.members, res.p_in);
  }
  BOOST_ASSERT_MSG(ranges::min(res.p_in) >= 0.0, "Implementation error: negative p emerged.");
  // Normalizes by sum
  auto sum = accumulate_sum(res.p_in);
  if (sum <= 0.0) {
    res.is_all_zero = true;
  } else {
    ranges::for_each(res.p_in, LAMBDA_1(_1 /= sum));
  }
  return res;
}

// The probability of spread success from gu.members[j0] to any member of gv
// is_seed: Whether the source vertex gu.members[j0] is a seed.
template <size_t N>
auto get_merged_p_simple(const CoarsenedEdgeDetails& edge, const CoarsenedVertexDetails& gu,
                         const CoarsenedVertexDetails& gv, size_t j0, bool is_seed) -> edge_probability_t {
  auto [p_cross_first, p_internal_first] =
      is_seed ? std::tie(edge.p_seed_cross, gu.p_seed_internal) : std::tie(edge.p_cross, gu.p_internal);

  if constexpr (N == 2) {
    // 2 : For each destination index i, at most 2 routes (see below)
    constexpr auto MAX_N_PATHS = 2 * CoarsenedVertexDetails::MAX_N_MEMBERS;
    auto paths = StaticVector<edge_probability_t, MAX_N_PATHS>{};
    auto j1 = get_next_j<2>(j0);
    for (auto i : gv.member_indices()) {
      paths.push_back(p_cross_first[j0][i]);                           // (1) j0 -> i
      paths.push_back(p_internal_first[j0][j1] * edge.p_cross[j1][i]); // (2) j0 -> j1 -> i
    }
    return at_least_1_probability_of_range(paths);
  } else if constexpr (N == 3) {
    // 5 : For each destination index i, at most 5 routes (see below)
    constexpr auto MAX_N_PATHS = 5 * CoarsenedVertexDetails::MAX_N_MEMBERS;
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
    return at_least_1_probability_of_range(paths);
  } else {
    static_assert(rfl::always_false_v<std::make_index_sequence<N>>, "Invalid group size.");
  }
}

using PreciseStateContainer = Array2D<            //
    edge_probability_t,                           //
    (1 << CoarsenedVertexDetails::MAX_N_MEMBERS), //
    (1 << CoarsenedVertexDetails::MAX_N_MEMBERS)  //
    >;

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
  constexpr auto MAX_N_PATHS = CoarsenedVertexDetails::MAX_N_MEMBERS * CoarsenedVertexDetails::MAX_N_MEMBERS;
  auto paths = StaticVector<edge_probability_t, MAX_N_PATHS>{};
  // Extracts all the direct paths in S -> Gv
  for (auto j : indices_in_bitset(N, S)) {
    for (auto i : gv.member_indices()) {
      paths.push_back(edge.p_cross[j][i]);
    }
  }
  auto p_success_in_this_turn = at_least_1_probability_of_range(paths);
  constexpr auto msg_pattern = "F[{1:#0{0}b}][{2:#0{0}b}]: Probability of this turn = {3:.4f}";
  // 2 : extra width of the prefix "0b"
  ELOGFMT(TRACE, msg_pattern, N + 2, T, S, p_success_in_this_turn);

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
      p_S_to_S_next *= at_least_1_probability_of_range( //
          indices_in_bitset(N, S) | TRANSFORM_VIEW(gu.p_internal[_1][dest]));
    }
    // Each direct path in S -> others must fail
    auto p_S_to_others = at_least_1_probability_of_range(
        views::cartesian_product(indices_in_bitset(N, S), indices_in_bitset(N, others)) | to_p_internal);
    // Message propagates to S_next in the next turn
    auto p_next = get_merged_p_precise_impl(edge, gu, gv, F, T | S_next, S_next);
    constexpr auto msg_pattern = "T = {1:#0{0}b}, S = {2:#0{0}b}, S_next = {3:#0{0}b}, others = {4:#0{0}b}, "
                                 "p_S_to_S_next = {5:.4f}, p_S_to_others = {6:.4f}, p_next = {7:.4f}";
    // 2 : Extra width of the prefix "0b"
    ELOGFMT(TRACE, msg_pattern, N + 2, T, S, S_next, others, p_S_to_S_next, p_S_to_others, p_next);
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

template <size_t N>
auto get_merged_p_from_single_source(const CoarsenedEdgeDetails& edge, const CoarsenedVertexDetails& gu,
                                     const CoarsenedVertexDetails& gv, size_t j0, EdgeWeightRule rule)
    -> edge_probability_t {
  if (rule == EdgeWeightRule::SEPARATE_SIMPLE || rule == EdgeWeightRule::MERGED_SIMPLE) {
    return get_merged_p_simple<N>(edge, gu, gv, j0, false);
  } else {
    return get_merged_p_precise(edge, gu, gv, make_bitset_from_indices({j0}));
  }
}

template <size_t N>
  requires(N == 2 || N == 3)
auto merge_coarsened_wim_graph_edge_common(CoarsenedEdgeDetails& dest, const CoarsenedEdgeDetails* inv_dest,
                                           const CoarsenedVertexDetails& gu, const CoarsenedVertexDetails& gv,
                                           const CoarseningParams& params) -> void {
  BOOST_ASSERT(gu.n_members() == N);
  // Part (1): merged.p
  auto [p_in, p_in_all_zero] = get_heuristic_p_in(dest, inv_dest, gu, gv, params.in_out_heuristic_rule);
  ELOGFMT(TRACE, "normalized in-heuristic of group {} = {::.4f}", gu.members, p_in);
  if (p_in_all_zero) {
    // Cannot serve as non-seed coarsened vertex since the group is isolated
    // (with no in-edge from other coarsened vertices)
    dest.merged.p = 0.0;
  } else if (params.edge_weight_rule != EdgeWeightRule::MERGED_PRECISE) {
    auto p_paths = std::array<edge_probability_t, N>{};
    for (auto j0 : range(N)) {
      p_paths[j0] = get_merged_p_from_single_source<N>(dest, gu, gv, j0, params.edge_weight_rule);
    }
    ELOGFMT(TRACE, "p_paths from group {} to group {} = {::.4f}", gu.members, gv.members, p_paths);

    auto to_combined = TRANSFORM_VIEW(p_in[_1] * p_paths[_1]);
    if (params.edge_weight_rule == EdgeWeightRule::MERGED_SIMPLE) {
      // Conditional probability, where VS is a virtual source
      // with out-edge VS -> gu.members[i] with p = p_in[i] for each i = 0 ... gu.n_members() - 1
      dest.merged.p = at_least_1_probability_of_range(range(N) | to_combined) / // P(VS -> gu.members -> gv.members)
                      at_least_1_probability_of_range(p_in);                    // P(VS -> gu.members)
    } else {
      BOOST_ASSERT(params.edge_weight_rule == EdgeWeightRule::SEPARATE_SIMPLE ||
                   params.edge_weight_rule == EdgeWeightRule::SEPARATE_PRECISE);
      // Weighted sum with p_in[i] as weight of gu.members[i]
      dest.merged.p = accumulate_sum(range(N) | to_combined);
    }
  } else {
    BOOST_ASSERT(params.edge_weight_rule == EdgeWeightRule::MERGED_PRECISE);
    constexpr auto U = (1zu << N) - 1;
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
      ELOGFMT(TRACE, msg_pattern, N + 2, S, p_enter_1, p_enter_0, p_pass); // 2 : extra width of the prefix "0b"
      sum += p_enter_1 * p_enter_0 * p_pass;
    }
    BOOST_ASSERT_MSG(sum <= 1.0, "Implementation error: probability sum is out of range [0, 1].");
    // Conditional probability: P(VS -> gu.members -> gv.members) / P(VS -> gu.members)
    dest.merged.p = sum / at_least_1_probability_of_range(p_in);
  }
  // Part (2): merged.p_boost, which only considers the best candidate inside the group gu
  auto estimated_gain_as_seed = std::array<edge_probability_t, N>{};
  for (auto j0 : range(N)) {
    estimated_gain_as_seed[j0] = get_merged_p_simple<N>(dest, gu, gv, j0, true);
  }
  auto best_at = ranges::max_element(estimated_gain_as_seed);
  ELOGFMT(TRACE, "Estimated gain as seed of group {}: {} (best vertex as seed: {})", //
          gu.members, estimated_gain_as_seed, gu.members[best_at - estimated_gain_as_seed.begin()]);
  dest.merged.p_seed = *best_at;
}

auto merge_coarsened_wim_graph_edge(CoarsenedEdgeDetails& dest, const CoarsenedEdgeDetails* inv_dest,
                                    const CoarsenedVertexDetails& gu, const CoarsenedVertexDetails& gv,
                                    const CoarseningParams& params) -> void {
  switch (gu.n_members()) {
  case 1: {
    auto [p_in, p_in_is_zero] = get_heuristic_p_in(dest, inv_dest, gu, gv, params.in_out_heuristic_rule);
    auto Mv = gv.n_members();
    dest.merged.p = p_in_is_zero ? 0.0_ep : at_least_1_probability_of_range(dest.p_cross[0] | views::take(Mv));
    dest.merged.p_seed = at_least_1_probability_of_range(dest.p_seed_cross[0] | views::take(Mv));
  } break;
  case 2:
    merge_coarsened_wim_graph_edge_common<2>(dest, inv_dest, gu, gv, params);
    break;
  case 3:
    merge_coarsened_wim_graph_edge_common<3>(dest, inv_dest, gu, gv, params);
    break;
  default:
    BOOST_ASSERT_MSG(false, "Unreachable branch.");
    std::abort();
  }
}

auto get_coarsened_wim_graph_edges(AdjacencyList<WIMEdge>& graph, const CoarseningDetails& details,
                                   const CoarseningParams& params) -> CoarseningDetails::EdgeDetailsMap {
  auto res = CoarseningDetails::EdgeDetailsMap{};
  // Step 1: builds p_cross relations
  for (auto [u, v, w] : graph::make_edge_range<0>(graph)) {
    auto [gu, gv] = details.to_group_id(u, v);
    if (gu == gv) {
      continue; // Skips the edges (u, v) whose two ends locate in the same group
    }
    auto it = res.find({gu, gv});
    if (it == res.end()) {
      auto [su, sv] = details.to_group_size(u, v);
      it = res.insert({{gu, gv}, make_initial_coarsened_edge_details(su, sv)}).first;
    }
    auto [ju, jv] = details.to_index_in_group(u, v);
    it->second.p_cross[ju][jv] = w.p;
    it->second.p_seed_cross[ju][jv] = w.p_seed;
  }
  // Step 2: gets merged p and p_boost
  for (auto& [key, dest] : res) {
    auto [gu, gv] = key;
    auto inv_dest = [&]() -> const CoarsenedEdgeDetails* {
      auto it = res.find({gv, gu});
      return it == res.end() ? nullptr : &it->second;
    }();
    merge_coarsened_wim_graph_edge(dest, inv_dest, details.groups[gu], details.groups[gv], params);
  }
  return res;
}

auto get_coarsened_wim_graph_vertex_weights(const CoarseningDetails& details,
                                            std::span<const vertex_weight_t> vertex_weights,
                                            const CoarseningParams& params) -> std::vector<vertex_weight_t> {
  auto coarsened_vertex_weights = std::vector<vertex_weight_t>(details.n_coarsened, 0.0);
  // vc = Coarsened vertex index in range [0, Nc)
  for (auto [vc, gc] : details.groups | views::enumerate) {
    auto h_in_sum = accumulate_sum(gc.heuristics_in);
    auto& dest = coarsened_vertex_weights[vc];

    if (params.vertex_weight_rule == VertexWeightRule::SUM) {
      dest = accumulate_sum(gc.members | views::transform(LAMBDA_1(vertex_weights[_1])));
    } else if (h_in_sum <= 0.0) {
      // Special case for AVERAGE and AVERAGE_BY_PATHS: Takes arithmetic average if no in-heuristic
      dest = accumulate_sum(gc.members | views::transform(LAMBDA_1(vertex_weights[_1]))) / gc.n_members();
    } else if (params.vertex_weight_rule == VertexWeightRule::AVERAGE) {
      // AVERAGE: Weighted average by normalized in-heuristic
      auto to_vertex_weight = LAMBDA_1(vertex_weights[gc.members[_1]] * gc.heuristics_in[_1] / h_in_sum);
      dest = accumulate_sum(range(gc.n_members()) | views::transform(to_vertex_weight));
    } else if (gc.members.size() == 1) {
      // AVERAGE_BY_PATHS for group of size 1
      dest = vertex_weights[gc.members[0]];
    } else if (gc.members.size() == 2) {
      // AVERAGE_BY_PATHS for group of size 2
      for (auto [i1, v1] : gc.members | views::enumerate) {
        // (1) v1; (2) v0 -> v1
        auto i0 = get_next_j<2>(i1);
        auto paths = {gc.heuristics_in[i1] / h_in_sum, gc.heuristics_in[i0] * gc.p_internal[i0][i1] / h_in_sum};
        dest += vertex_weights[v1] * at_least_1_probability_of_range(paths);
      }
    } else {
      // AVERAGE_BY_PATHS for group of size 3
      BOOST_ASSERT(gc.members.size() == 3);
      for (auto [i2, v2] : gc.members | views::enumerate) {
        auto [i0, i1] = get_next_j<3>(i2);
        auto prob = at_least_1_probability(
            gc.heuristics_in[i2] / h_in_sum,                         // (1)   v2
            gc.heuristics_in[i0] * gc.p_internal[i0][i2] / h_in_sum, // (2,3) v0 -> v2, or v0 -> v1 -> v2
            gc.heuristics_in[i0] * gc.p_internal[i0][i1] * gc.p_internal[i1][i2] / h_in_sum,
            gc.heuristics_in[i1] * gc.p_internal[i1][i2] / h_in_sum, // (4,5) v1 -> v2, or v1 -> v0 -> v2
            gc.heuristics_in[i1] * gc.p_internal[i1][i0] * gc.p_internal[i0][i2] / h_in_sum);
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
  ELOGFMT(DEBUG, "Parameters: {:4}", params);
  auto n = graph::num_vertices(graph);
  // d_in: unit weight by default
  auto details = CoarseningDetails{.n = n,
                                   .n_coarsened = n_groups,
                                   .group_id = std::vector(group_id.begin(), group_id.end()),
                                   .index_in_group = std::vector<vertex_id_t>(n),
                                   .groups = std::vector<CoarsenedVertexDetails>(n_groups)};
  ELOGFMT(DEBUG, "Starts coarsening: n = {}, n_groups = {}", details.n, details.n_coarsened);
  // Step 1: groups[].members & groups[].vertex_weights & groups[].heuristics_in & index_in_group
  // Vertex v is placed in group g
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
  ELOG_DEBUG << "Done coarsening step 1: assigning group members.";
  // Step 2: groups[].p_internal
  for (auto [u, v, w] : graph::make_edge_range<0>(graph)) {
    auto g = group_id[u];
    if (g == group_id[v]) {
      auto [ju, jv] = details.to_index_in_group(u, v);
      details.groups[g].p_internal[ju][jv] = w.p;
      details.groups[g].p_seed_internal[ju][jv] = w.p_seed;
    }
  }
  ELOG_DEBUG << "Done coarsening step 2: assigning p_internal of groups.";
  // Step 3: Generates coarsened edges
  details.edges = get_coarsened_wim_graph_edges(graph, details, params);
  ELOG_DEBUG << "Done coarsening step 3: coarsening edges.";
  // Step 4: Generates coarsened vertex weights
  auto coarsened_vertex_weights = get_coarsened_wim_graph_vertex_weights(details, vertex_weights, params);
  ELOG_DEBUG << "Done coarsening step 4: coarsening vertex weights.";
  // Finally, builds the result from all the components
  auto coarsened_edge_list = DirectedEdgeList<WIMEdge>{};
  coarsened_edge_list.open_for_push_back();
  for (const auto& [p, w] : details.edges) {
    coarsened_edge_list.push_back(get<0>(p), get<1>(p), w.merged);
  }
  coarsened_edge_list.close_for_push_back();

  auto res_details = [&]() {
    if constexpr (std::is_same_v<ResultType, CoarsenGraphDetailedResult>) {
      // Returns the details directly
      return std::move(details);
    } else if constexpr (std::is_same_v<ResultType, CoarsenGraphBriefResult>) {
      // Brief information
      auto brief = CoarseningBrief{.n = n, .n_coarsened = n_groups, .groups = {}};
      brief.groups.reserve(n_groups);
      for (const auto& u : details.groups) {
        auto best_seed_index = select_best_seed_in_group(u).index_in_group;
        brief.groups.push_back({.members = u.members, .best_seed_index = best_seed_index});
      }
      return brief;
    } else {
      static_assert(rfl::always_false_v<ResultType>, "Invalid result type.");
    }
  }();
  return ResultType{.coarsened_graph = AdjacencyList<WIMEdge>{coarsened_edge_list},
                    .coarsened_inv_graph = InvAdjacencyList<WIMEdge>{coarsened_edge_list},
                    .coarsened_vertex_weights = std::move(coarsened_vertex_weights),
                    .details = std::move(res_details)};
}

template <class ResultType>
auto coarsen_wim_graph_w_framework(const AdjacencyList<WIMEdge>& graph, const InvAdjacencyList<WIMEdge>& inv_graph,
                                   std::span<const vertex_weight_t> vertex_weights,
                                   const CoarseningParams& params) noexcept -> ResultType {
  // Step 1:
  ELOGFMT(DEBUG, "|V|, |E| of graph = {}; |V|, |E| of inv_graph = {}; Size of vertex_weights = {}", //
          graph_n_m(graph), graph_n_m(inv_graph), vertex_weights.size());
  auto bidir_graph = merge_wim_edge_to_undirected(graph, params);
  ELOGFMT(DEBUG, "|V|, |E| of bidir_graph = {}", graph_n_m(bidir_graph));
  // Step 2
  auto [n_groups, group_id] = mongoose_match(bidir_graph, params);
  ELOGFMT(DEBUG, "n_groups = {}; Size of group_id = {}", n_groups, group_id.size());
  // Step 3
  return coarsen_wim_graph_w_impl<ResultType>(const_cast<AdjacencyList<WIMEdge>&>(graph),        //
                                              const_cast<InvAdjacencyList<WIMEdge>&>(inv_graph), //
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
      auto best_index = [&]() {
        if constexpr (std::is_same_v<CoarsenResultType, CoarsenGraphBriefResult>) {
          return group.best_seed_index;
        } else if constexpr (std::is_same_v<CoarsenResultType, CoarsenGraphDetailedResult>) {
          return select_best_seed_in_group(group).index_in_group;
        } else {
          static_assert(rfl::always_false_v<CoarsenResultType>, "Invalid coarsening result type.");
        }
      }();
      expanded_seeds.push_back(group.members[best_index]);
    }
  };
  ELOGFMT(DEBUG, "coarsened_seeds = {}", coarsened_seeds);

  if (params.seed_expanding_rule == SeedExpandingRule::LOCAL) {
    // S_LOCAL: Only considers information inside current group
    make_expanded_seeds_as_s_local();
  } else if (params.seed_expanding_rule == SeedExpandingRule::SIMULATIVE) {
    // S_SIMULATIVE: Picks one-by-one separately by simulation result of single expanded seed vertex
    for (auto sc : coarsened_seeds) {
      auto& group = coarsening_result.details.groups[sc];
      auto estimated_gain = StaticVector<double, CoarsenedVertexDetails::MAX_N_MEMBERS>{};
      for (auto v : group.members) {
        auto sim_res = wim_simulate_w(graph, vertex_weights, {n, {v}}, **params.simulation_try_count);
        RFL_RETURN_ON_ERROR(sim_res);
        estimated_gain.push_back(*sim_res);
      }
      auto best_expanded_seed = group.members[ranges::max_element(estimated_gain) - estimated_gain.begin()];
      constexpr auto msg_pattern = "Expanding seeds: Estimated gain of group {} "
                                   "(coarsened to vertex #{}) = {::.4f} (best member: {})";
      ELOGFMT(DEBUG, msg_pattern, group.members, sc, estimated_gain, best_expanded_seed);
      expanded_seeds.push_back(best_expanded_seed);
    }
  } else {
    // S_ITERATIVE: Result of S_SEPARATE as initial solution
    BOOST_ASSERT(params.seed_expanding_rule == SeedExpandingRule::ITERATIVE);
    make_expanded_seeds_as_s_local();
    auto sim_res = wim_simulate_w_s(graph, vertex_weights, expanded_seeds, **params.simulation_try_count);
    if (!sim_res) {
      return *sim_res.error(); // Handles unexpected failure
    }
    ELOGFMT(DEBUG, "Initial expanded seed set: {}\n\tWith simulation result = {:.4f}", expanded_seeds, *sim_res);
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
      ELOGFMT(DEBUG, "Iteration #{}: Attempting to change {} to {} (in group {} coarsened to vertex #{})",
              iteration_index, prev, next, group.members, coarsened_seeds[si]);
      expanded_seeds[si] = next;
      auto next_sim_res = wim_simulate_w_s(graph, vertex_weights, expanded_seeds, **params.simulation_try_count);
      if (!next_sim_res) {
        return *next_sim_res.error();
      }
      ELOGFMT(DEBUG, "Iteration #{}: Candidate {} vs. {}: Simulation result: {:.4f} vs. {:.4f}", iteration_index, next,
              prev, *next_sim_res, *sim_res);
      if (*next_sim_res > *sim_res) {
        sim_res = next_sim_res; // Changes prev to next
      } else {
        expanded_seeds[si] = prev; // Otherwise, restores current expanded seed
      }
      ELOGFMT(DEBUG, "Expanded seed set after iteration #{}: {}", iteration_index, expanded_seeds);
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
  ELOGFMT(DEBUG, "Starts Mongoose matching algorithm. |V| = {}", n);

  auto match = std::vector<vertex_id_t>(n, -1);
  auto adopt = std::vector<vertex_id_t>(n, -1);
  auto make_match = [&](vertex_id_t u, vertex_id_t v) {
    match[u] = v;
    match[v] = u;
  };

  // Step 1: Neighbor matching
  ranges::shuffle(seq, rand_engine);
  if (params.seed_merging_rule == SeedMergingRule::SINGLE) {
    for (auto s : seeds) {
      BOOST_ASSERT_MSG(s >= 0 && s < n, "Seed vertex out of range [0, n).");
      match[s] = s; // Marks the seeds as single-matched
    }
  }
  for (auto v : seq) {
    if (match[v] != -1) {
      continue; // Skips the matched vertices (including seeds)
    }
    auto best_p = initial_p(params.neighbor_match_rule);
    auto best_u = -1_vid;
    for (auto [u, p] : graph[v]) {
      // Selects an unmatched neighbor with the best edge weight
      if (match[u] == -1 && is_better_p(best_p, p, params.neighbor_match_rule)) {
        best_p = p;
        best_u = u;
      }
    }
    if (best_u == -1) {
      match[v] = v;
    } else {
      make_match(v, best_u);
    }
  }
  ELOG_TRACE << "Step 1 done: Neighbor matching.\n" << DUMP_INDEX_ARRAY(match);

  // Step 2: Brotherly matching & Adoptive matching
  ranges::shuffle(seq, rand_engine);
  if (params.seed_merging_rule == SeedMergingRule::SINGLE) {
    for (auto s : seeds) {
      BOOST_ASSERT_MSG(s >= 0 && s < n, "Seed vertex out of range [0, n).");
      match[s] = -1; // Uses -1 to mark seeds
    }
  }
  for (auto v : seq) {
    if (match[v] != v) {
      continue; // Skips the matched vertices or seeds
    }
    auto v_pivot = [&] {
      auto best_p = initial_p(params.neighbor_match_rule);
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
    ELOGFMT(TRACE, "v = {}: pivot = {}", v, v_pivot);
    if (v_pivot == -1) {
      ELOGFMT(TRACE, "Found isolated vertex {} during brotherly matching.", v);
      continue; // Special case when v is an isloated vertex, or all the neighbors are seeds
    }
    // Step 2.1 for each vertex v: Neighbor matching
    auto u_last = -1_vid;
    for (auto u : graph[v_pivot] | views::keys | views::filter(LAMBDA_1(match[_1] == _1))) {
      if (u_last != -1) {
        ELOGFMT(TRACE, "Brotherly matching around pivot #{}: {} -- {}", v_pivot, u_last, u);
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
    ELOGFMT(TRACE, "Begins adoptive matching for #{}, with pivot = {}, match[pivot] = {}", u_last, v_pivot, x);
    // v1 -> v2 : v2 adopts v1;  v1 -- v2 : (v1, v2) forms a match
    if (adopt[v_pivot] != -1) {
      constexpr auto msg_pattern_1 = "Adoptive matching (1): {0} -> {1} -- {2}  ==>  {0} -- {1}, {2} -- {3}";
      ELOGFMT(TRACE, msg_pattern_1, adopt[v_pivot], v_pivot, x, u_last);
      make_match(u_last, x);
      make_match(adopt[v_pivot], v_pivot);
      adopt[v_pivot] = -1;
    } else if (adopt[x] != -1) {
      constexpr auto msg_pattern_2 = "Adoptive matching (2): {0} -> {1} -- {2}  ==>  {0} -- {1}, {2} -- {3}";
      ELOGFMT(TRACE, msg_pattern_2, adopt[x], x, v_pivot, u_last);
      make_match(u_last, v_pivot);
      make_match(adopt[x], x);
      adopt[x] = -1;
    } else {
      constexpr auto msg_pattern_3 = "Adoptive matching (3): {0} -- {1}  ==>  {0} -- {1} <- {2}";
      ELOGFMT(TRACE, msg_pattern_3, x, v_pivot, u_last);
      match[u_last] = v_pivot;
      adopt[v_pivot] = u_last;
    }
  }
  ELOG_TRACE << "Step 2 done: Brotherly matching & Adoptive matching.\n"
             << DUMP_INDEX_ARRAY(match) << '\n'
             << DUMP_INDEX_ARRAY(adopt);

  // Step 3: Groups by matching result
  auto res = MongooseMatchResult{.n_groups = 0, .group_id = std::vector<vertex_id_t>(n, -1)};
  for (auto v : range(n)) {
    if (res.group_id[v] != -1) {
      continue; // Ignores the already grouped vertices
    }
    if (match[v] == -1 || match[v] == v) {
      // (1) v is isolated, or v is a seed vertex
      ELOGFMT(TRACE, "Isolated vertex found during grouping: {}", v);
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

auto further_coarsen_wim_graph(const CoarsenGraphBriefResult& last_result, const CoarseningParams& params) noexcept
    -> CoarsenGraphBriefResult {
  return coarsen_wim_graph_w(last_result.coarsened_graph, last_result.coarsened_inv_graph,
                             last_result.coarsened_vertex_weights, params);
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

auto further_coarsen_wim_graph_d(const CoarsenGraphBriefResult& last_result, const CoarseningParams& params) noexcept
    -> CoarsenGraphDetailedResult {
  return coarsen_wim_graph_d_w(last_result.coarsened_graph, last_result.coarsened_inv_graph,
                               last_result.coarsened_vertex_weights, params);
}

// ---- Step 4 ----

auto select_best_seed_in_group(const CoarsenedVertexDetails& v) noexcept -> SelectBestSeedResult {
  auto m = v.members.size();
  if (m == 1) {
    return {.index_in_group = 0};
  }
  BOOST_ASSERT(m <= 3);
  auto to_estimated_gain = [&](size_t j0) {
    auto res = edge_probability_t{};
    if (v.n_members() == 2) {
      auto j1 = get_next_j<2>(j0);
      res = v.heuristics_out_seed[j0] + v.p_seed_internal[j0][j1] * v.heuristics_out[j1];
    } else {
      BOOST_ASSERT(v.n_members() == 3);
      auto [j1, j2] = get_next_j<3>(j0);
      auto p1 = at_least_1_probability(                      //
          v.p_seed_internal[j0][j1],                         // j0 -> j1
          v.p_seed_internal[j0][j2] * v.p_internal[j2][j1]); // j0 -> j2 -> j1
      auto p2 = at_least_1_probability(                      //
          v.p_seed_internal[j0][j2],                         // j0 -> j2
          v.p_seed_internal[j0][j1] * v.p_internal[j1][j2]); // j0 -> j1 -> j2
      res = v.heuristics_out_seed[j0] + p1 * v.heuristics_out[j1] + p2 * v.heuristics_out[j2];
    }
    return res;
  };
  auto estimated_gain = [&]() {
    auto view = range(m) | views::transform(to_estimated_gain);
    return std::vector(view.begin(), view.end());
  }();
  ELOGFMT(TRACE, "Best seed selection: estimated gain for group {} = {::.4f}", //
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
  return expand_wim_seed_vertices_w( //
      last_result.coarsened_graph, last_result.coarsened_vertex_weights, cur_result, coarsened_seeds, params);
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
  return expand_wim_seed_vertices_d_w( //
      last_result.coarsened_graph, last_result.coarsened_vertex_weights, cur_result, coarsened_seeds, params);
}

// ---- Dump functions ----

namespace {
template <ranges::forward_range Range, class ToStringFn>
  requires(std::is_invocable_r_v<std::string, ToStringFn, ranges::range_value_t<Range>>)
auto dump_array_as_braced_list(Range&& values, ToStringFn&& to_str_fn, int indent, int level) -> std::string {
  constexpr auto N_VALUES_PER_ROW = 10;
  if (indent <= 0) {
    auto joined = values | views::transform(to_str_fn) | views::join_with(", "s) | views::common;
    return '{' + std::string(joined.begin(), joined.end()) + '}';
  }
  auto res = "{"s;
  for (auto [i, x] : values | views::enumerate) {
    if (i != 0) {
      res += ", ";
    }
    if (i % N_VALUES_PER_ROW == 0) {
      res += '\n' + std::string(indent * (level + 1), ' ');
    }
    res += to_str_fn(x);
  }
  res += '\n' + std::string(indent * level, ' ') + '}';
  return res;
}

template <ranges::forward_range Range>
auto dump_array_as_braced_list(Range&& values, int width, int indent, int level) -> std::string {
  auto to_str_fn = LAMBDA_1(fmt::format("{0:{1}}", _1, width));
  return dump_array_as_braced_list(std::forward<Range>(values), to_str_fn, indent, level);
}

template <forward_range_of<std::string> ComponentRange>
auto merge_dumped_components(ComponentRange&& components, int indent, int level) -> std::string {
  if (indent <= 0) {
    auto joined = components | views::join_with(", "s) | views::common;
    return '{' + std::string(joined.begin(), joined.end()) + '}';
  }
  auto res = "{"s;
  for (const auto& [i, c] : components | views::enumerate) {
    if (i != 0) {
      res += ',';
    }
    res += '\n' + std::string(indent * (level + 1), ' ');
    res += c;
  }
  res += '\n' + std::string(indent * level, ' ') + '}';
  return res;
}

template <is_edge_property E, int IsInv>
auto dump_graph_generic(const graph::adjacency<IsInv, E>& graph, int indent, int level) -> std::string {
  auto edge_range = graph::make_edge_range<0>(remove_const(graph));
  auto components = make_reserved_vector<std::string>(graph.num_edges());
  for (auto [u, v, w] : edge_range) {
    components.push_back(fmt::format("({}, {}): {}", u, v, w));
  }
  return merge_dumped_components(components, indent, level);
}

template <class CoarseningDetailsType>
auto dump_impl(const CoarsenGraphResult<CoarseningDetailsType>& result_obj, int indent, int level) -> std::string {
  constexpr auto DECIMAL_DIGITS = 4;
  auto weights_width = // 1 : One position for the decimal point '.'
      static_cast<int>(std::log10(ranges::max(result_obj.coarsened_vertex_weights))) + DECIMAL_DIGITS + 1;
  auto weight_to_str = LAMBDA_1(fmt::format("{0:{1}.{2}f}", _1, weights_width, DECIMAL_DIGITS));
  auto weights_str = dump_array_as_braced_list(result_obj.coarsened_vertex_weights, weight_to_str, indent, level + 1);
  auto components = {
      fmt::format(".coarsened_graph = {}", dump_graph_generic(result_obj.coarsened_graph, indent, level + 1)),
      fmt::format(".coarsened_inv_graph = {}", dump_graph_generic(result_obj.coarsened_inv_graph, indent, level + 1)),
      fmt::format(".coarsened_vertex_weights = {}", weights_str),
      fmt::format(".details = {}", result_obj.details.dump(indent, level + 1))};
  return merge_dumped_components(components, indent, level);
}
} // namespace

auto CoarsenedVertexBrief::dump(int indent, int level) const noexcept -> std::string {
  indent = std::max(indent, 0);
  level = std::max(level, 0);
  auto components = {fmt::format(".members = {}", members), fmt::format(".best_seed_index = {}", best_seed_index)};
  return merge_dumped_components(components, indent, level);
}

auto CoarsenedVertexDetails::dump(int indent, int level) const noexcept -> std::string {
  indent = std::max(indent, 0);
  level = std::max(level, 0);
  auto dump_p_internal = [&](const PInternalContainer& p_internal) {
    auto m = n_members();
    auto rows = range(m) | views::transform(LAMBDA_1(fmt::format("{::.4f}", p_internal[_1] | views::take(m))));
    if (indent <= 0) {
      return fmt::format("{}", rows);
    }
    auto res = "["s;
    for (auto i : range(m)) {
      res += '\n' + std::string(indent * (level + 2), ' ');
      res += rows[i];
    }
    res += '\n' + std::string(indent * (level + 1), ' ') + ']';
    return res;
  };
  auto components = {fmt::format(".members = {}", members),
                     fmt::format(".vertex_weights = {::.4f}", vertex_weights),
                     fmt::format(".heuristics_in = {::.4f}", heuristics_in),
                     fmt::format(".heuristics_out = {::.4f}", heuristics_out),
                     fmt::format(".heuristics_out_seed = {::.4f}", heuristics_out_seed),
                     fmt::format(".p_internal = {}", dump_p_internal(p_internal)),
                     fmt::format(".p_seed_internal = {}", dump_p_internal(p_seed_internal))};
  return merge_dumped_components(components, indent, level);
}

auto CoarsenedEdgeDetails::dump(int indent, int level) const noexcept -> std::string {
  indent = std::max(indent, 0);
  level = std::max(level, 0);
  auto dump_p_cross = [&](const PCrossContainer& p_cross) {
    auto n_rows = n_members_left;
    auto n_cols = n_members_right;
    auto rows = range(n_rows) | views::transform(LAMBDA_1(fmt::format("{::.4f}", p_cross[_1] | views::take(n_cols))));
    if (indent <= 0) {
      return fmt::format("{}", rows);
    }
    auto res = "["s;
    for (auto i : range(n_rows)) {
      res += '\n' + std::string(indent * (level + 2), ' ');
      res += rows[i];
    }
    res += '\n' + std::string(indent * (level + 1), ' ') + ']';
    return res;
  };
  auto components = {fmt::format(".p_cross = {}", dump_p_cross(p_cross)),           //
                     fmt::format(".p_seed_cross = {}", dump_p_cross(p_seed_cross)), //
                     fmt::format(".merged = {}", merged)};
  return merge_dumped_components(components, indent, level);
}

auto CoarseningBrief::dump(int indent, int level) const noexcept -> std::string {
  indent = std::max(indent, 0);
  level = std::max(level, 0);
  auto groups_str = [&]() {
    auto to_group_str = LAMBDA_1(fmt::format("[{}] = {}", _1, ::dump(groups[_1])));
    return merge_dumped_components(range(n_coarsened) | views::transform(to_group_str), indent, level + 1);
  }();
  auto components = {fmt::format(".n = {}", n), fmt::format(".n_coarsened = {}", n_coarsened),
                     fmt::format(".groups = {}", groups_str)};
  return merge_dumped_components(components, indent, level);
}

auto CoarseningDetails::dump(int indent, int level) const noexcept -> std::string {
  indent = std::max(indent, 0);
  level = std::max(level, 0);
  auto groups_str = [&]() {
    auto to_group_str = LAMBDA_1(fmt::format("[{}] = {}", _1, groups[_1].dump(indent, level + 2)));
    return merge_dumped_components(range(n_coarsened) | views::transform(to_group_str), indent, level + 1);
  }();
  auto edges_str = [&]() {
    auto to_edge_str = [&](const EdgeDetailsMap::value_type& pair) {
      const auto& [v_pair, e] = pair;
      return fmt::format("{} = {}", v_pair, e.dump(indent, level + 2));
    };
    return merge_dumped_components(edges | views::transform(to_edge_str), indent, level + 1);
  }();
  auto group_id_width = static_cast<int>(std::log10(n_coarsened)) + 1;
  auto components = {
      fmt::format(".n = {}", n),
      fmt::format(".n_coarsened = {}", n_coarsened),
      fmt::format(".group_id = {}", dump_array_as_braced_list(group_id, group_id_width, indent, level + 1)),
      fmt::format(".index_in_group = {}", dump_array_as_braced_list(index_in_group, 1, indent, level + 1)),
      fmt::format(".groups = {}", groups_str),
      fmt::format(".edges = {}", edges_str)};
  return merge_dumped_components(components, indent, level);
}

template <>
auto CoarsenGraphResult<CoarseningBrief>::dump(int indent, int level) const noexcept -> std::string {
  return dump_impl(*this, indent, level);
}

template <>
auto CoarsenGraphResult<CoarseningDetails>::dump(int indent, int level) const noexcept -> std::string {
  return dump_impl(*this, indent, level);
}
