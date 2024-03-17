#include "coarsening.h"
#include "dump.h"
#include "utils/easylog.h"
#include "wim.h"
#include <fmt/ranges.h>
#include <magic_enum.hpp>
#include <nwgraph/adaptors/edge_range.hpp>

namespace {
constexpr auto MAX_N_MEMBERS = CoarsenedVertexDetailsBase::MAX_N_MEMBERS;
constexpr auto P_NOT_ASSIGNED = -1.0_ep;

using VertexPair = std::pair<vertex_id_t, vertex_id_t>;
using HeuristicsContainer = CoarsenedVertexDetailsBase::HeuristicsContainer;
using PCrossContainer = CoarsenedEdgeDetailsBase::PCrossContainer;
using PInternalContainer = CoarsenedVertexDetailsBase::PInternalContainer;

struct StepTimer {
  std::string task_name;
  nw::util::seconds_timer timer;
  int step_count = 0;
  double total_time_used;

  explicit StepTimer(std::string task_name)
      : task_name(std::move(task_name)), timer(), step_count(0), total_time_used(0.0) {}

  template <class Func, class... Args>
    requires(std::is_invocable_v<Func, Args...>)
  auto do_step(std::string_view step_description, Func&& func, Args&&... args) {
    constexpr auto RETURNS_VOID = std::is_void_v<std::invoke_result_t<Func, Args...>>;

    timer.start();
    auto invoke_result = [&]() {
      if constexpr (RETURNS_VOID) {
        std::invoke(func, std::forward<Args>(args)...);
        return rfl::Nothing{};
      } else {
        return std::invoke(func, std::forward<Args>(args)...);
      }
    }();
    timer.stop();

    step_count += 1;
    total_time_used += timer.elapsed();
    constexpr auto msg_pattern = "Done step {} of task '{}': {}."
                                 "\n\tTime usage: {:.3f} in current step, {:.3f} in total.";
    ELOGFMT(INFO, msg_pattern, step_count, task_name, step_description, timer.elapsed(), total_time_used);
    if constexpr (!RETURNS_VOID) {
      return invoke_result;
    }
  }
};

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

template <auto WhichP, is_edge_property E, int IsInv,
          same_as_either<std::nullptr_t, VertexSet> SkipSet = std::nullptr_t>
  requires(member_object_pointer_from_to<decltype(WhichP), E, edge_probability_t>)
auto get_heuristic_generic(const graph::adjacency<IsInv, E>& graph, std::span<const vertex_id_t> group_id,
                           std::span<const vertex_weight_t> vertex_weights, vertex_id_t v, InOutHeuristicRule rule,
                           const SkipSet& skip_set = nullptr) -> edge_probability_t {
  auto res = 0.0_ep;
  for (auto [u, w] : graph[v]) {
    if (group_id[u] == group_id[v]) {
      continue;
    }
    if constexpr (std::is_same_v<SkipSet, VertexSet>) {
      if (skip_set.contains(u)) {
        continue; // Customized skipping rule
      }
    }
    if (rule == InOutHeuristicRule::P) {
      res += w.*WhichP;
    } else if (rule == InOutHeuristicRule::W) {
      res += w.*WhichP * vertex_weights[u];
    } else {
      BOOST_ASSERT(rule == InOutHeuristicRule::COUNT || rule == InOutHeuristicRule::UNIT);
      res += 1.0_ep; // Counts only
    }
  }
  return res;
}

#define DEFINE_GET_HEURISTIC_FUNCTION(suffix, IsInv, p_member)                                                 \
  template <is_edge_property E>                                                                                \
  auto get_heuristics_##suffix(const graph::adjacency<IsInv, E>& graph, std::span<const vertex_id_t> group_id, \
                               std::span<const vertex_weight_t> vertex_weights, vertex_id_t v,                 \
                               InOutHeuristicRule rule) -> edge_probability_t {                                \
    return get_heuristic_generic<&E::p_member>(graph, group_id, vertex_weights, v, rule);                      \
  }                                                                                                            \
  template <is_edge_property E>                                                                                \
  auto get_heuristics_##suffix(const graph::adjacency<IsInv, E>& graph, std::span<const vertex_id_t> group_id, \
                               std::span<const vertex_weight_t> vertex_weights, const VertexSet& skip_set,     \
                               vertex_id_t v, InOutHeuristicRule rule) -> edge_probability_t {                 \
    return get_heuristic_generic<&E::p_member>(graph, group_id, vertex_weights, v, rule, skip_set);            \
  }

DEFINE_GET_HEURISTIC_FUNCTION(in, 1, p)
DEFINE_GET_HEURISTIC_FUNCTION(in_boost, 1, p_boost)
DEFINE_GET_HEURISTIC_FUNCTION(out, 0, p)
DEFINE_GET_HEURISTIC_FUNCTION(out_seed, 0, p_seed)

#undef DEFINE_GET_HEURISTIC_FUNCTION

struct NormalizedHeuristics {
  enum PolicyWhenAllZero { KEEPS_WHEN_ALL_ZERO, TO_UNIFORM_WHEN_ALL_ZERO };

  HeuristicsContainer p;
  bool is_all_zero;

  NormalizedHeuristics(HeuristicsContainer p_raw, InOutHeuristicRule rule, PolicyWhenAllZero policy)
      : p(std::move(p_raw)) {
    constexpr auto EPS = 1.0e-8_ep;
    BOOST_ASSERT_MSG(p.size() > 0, "Empty group disallowed.");
    // Ensures that each value falls in range [0, +inf) in case of floating point error
    for (auto& h : p) {
      BOOST_ASSERT(!std::isinf(h) && !std::isnan(h));
      h = std::max(h, 0.0_ep);
    }
    // For UNIT rule, normalizes to {0, 1} first.
    if (rule == InOutHeuristicRule::UNIT) {
      ranges::for_each(p, LAMBDA_1(_1 = (_1 < EPS) ? 0.0_ep : 1.0_ep));
    }
    // Normalizes by sum such that accumulate_sum(p_in) == 1.0
    auto sum = accumulate_sum(p);
    if (sum <= 0.0) {
      if (policy == KEEPS_WHEN_ALL_ZERO) {
        is_all_zero = true;
      } else {
        BOOST_ASSERT(policy == TO_UNIFORM_WHEN_ALL_ZERO);
        ranges::fill(p, 1.0_ep / p.size()); // Sets all the heuristic values to 1 / M
        is_all_zero = false;
      }
    } else {
      ranges::for_each(p, LAMBDA_1(_1 /= sum));
      is_all_zero = false;
    }
  }
};

#define DEFINE_GET_NORMALIZED_HEURISTIC_FUNCTION(suffix, heuristic_member)                                         \
  template <is_coarsened_vertex_details DetailsType>                                                               \
    requires(requires(DetailsType x) {                                                                             \
      { x.heuristic_member } -> std::same_as<HeuristicsContainer&>;                                                \
    }) /* Using C++20 constraint to check member presence */                                                       \
  auto get_normalized_heuristic_##suffix(const DetailsType& gu, InOutHeuristicRule rule,                           \
                                         NormalizedHeuristics::PolicyWhenAllZero policy) -> NormalizedHeuristics { \
    return NormalizedHeuristics(gu.heuristic_member, rule, policy);                                                \
  }

DEFINE_GET_NORMALIZED_HEURISTIC_FUNCTION(p_in, heuristics_in)
DEFINE_GET_NORMALIZED_HEURISTIC_FUNCTION(p_in_boost, heuristics_in_boost)
DEFINE_GET_NORMALIZED_HEURISTIC_FUNCTION(p_out, heuristics_out)
DEFINE_GET_NORMALIZED_HEURISTIC_FUNCTION(p_out_seed, heuristics_out_seed)

#undef DEFINE_GET_NORMALIZED_HEURISTIC_FUNCTION

template <is_coarsened_edge_details EdgeDetailsType, is_coarsened_vertex_details VertexDetailsType>
auto get_normalized_heuristic_p_in(const EdgeDetailsType& edge, const EdgeDetailsType* inv_edge,
                                   const VertexDetailsType& gu, const VertexDetailsType& gv, InOutHeuristicRule rule,
                                   NormalizedHeuristics::PolicyWhenAllZero policy) -> NormalizedHeuristics {
  constexpr auto EPS = 1.0e-8_ep;
  auto p_temp = gu.heuristics_in;
  if (inv_edge != nullptr) {
    // Removes the directed edge gu.members[j] <-- gv.members[i]
    for (auto [j, i] : views::cartesian_product(gu.member_indices(), gv.member_indices())) {
      if (rule == InOutHeuristicRule::P) {
        p_temp[j] -= inv_edge->p_cross[i][j];
      } else if (rule == InOutHeuristicRule::W) {
        p_temp[j] -= inv_edge->p_cross[i][j] * gv.vertex_weights[i];
      } else {
        BOOST_ASSERT(rule == InOutHeuristicRule::UNIT || rule == InOutHeuristicRule::COUNT);
        p_temp[j] -= (inv_edge->p_cross[i][j] < EPS) ? 0.0_ep : 1.0_ep;
      }
    }
    constexpr auto msg_pattern = "In-heuristics of group {} after removing back edges from {} = {::.4f}";
    MYLOG_FMT_TRACE(msg_pattern, gu.members, gv.members, p_temp);
  }
  return NormalizedHeuristics(p_temp, rule, policy);
}

// The probability of spread success from Gu.members[j0] to Gv.members[i]
auto get_merged_p_simple_to_with_first(const PCrossContainer& edge_p_cross, const PCrossContainer& edge_p_cross_first,
                                       const PInternalContainer& gU_p_internal,
                                       const PInternalContainer& gU_p_internal_first, //
                                       size_t gU_n_members, size_t j0, size_t i) -> edge_probability_t {
  if (gU_n_members == 1) {
    return edge_p_cross_first[j0][i];
  } else if (gU_n_members == 2) {
    auto j1 = get_next_j<2>(j0);
    return at_least_1_probability(                        //
        edge_p_cross_first[j0][i],                        // j0 -> i
        gU_p_internal_first[j0][j1] * edge_p_cross[j1][i] // j0 -> j1 -> i
    );
  }
  BOOST_ASSERT_MSG(gU_n_members == 3, "Group size other than 1, 2, 3 is not supported.");
  auto [j1, j2] = get_next_j<3>(j0);
  return at_least_1_probability(                                                 //
      edge_p_cross_first[j0][i],                                                 // j0 -> i
      gU_p_internal_first[j0][j1] * edge_p_cross[j1][i],                         // j0 -> j1 -> i
      gU_p_internal_first[j0][j2] * edge_p_cross[j2][i],                         // j0 -> j2 -> i
      gU_p_internal_first[j0][j1] * gU_p_internal[j1][j2] * edge_p_cross[j2][i], // j0 -> j1 -> j2 -> i
      gU_p_internal_first[j0][j2] * gU_p_internal[j2][j1] * edge_p_cross[j1][i]  // j0 -> j2 -> j1 -> i
  );
}

auto get_merged_p_simple_to(const PCrossContainer& edge_p_cross, const PInternalContainer& gU_p_internal,
                            size_t gU_n_members, size_t j0, size_t i) -> edge_probability_t {
  return get_merged_p_simple_to_with_first( //
      edge_p_cross, edge_p_cross, gU_p_internal, gU_p_internal, gU_n_members, j0, i);
}

// The probability of spread success from gu.members[j0] to any member of gv
auto get_merged_p_simple_with_first(const PCrossContainer& edge_p_cross, const PCrossContainer& edge_p_cross_first,
                                    const PInternalContainer& gU_p_internal,
                                    const PInternalContainer& gU_p_internal_first, size_t gU_n_members,
                                    size_t gV_n_members, size_t j0) -> edge_probability_t {
  auto to_p_to_target = TRANSFORM_VIEW(get_merged_p_simple_to_with_first( //
      edge_p_cross, edge_p_cross_first, gU_p_internal, gU_p_internal_first, gU_n_members, j0, _1));
  return at_least_1_probability_r(range(gV_n_members) | to_p_to_target);
}

// The probability of spread success from gu.members[j0] to any member of gv
auto get_merged_p_simple(const PCrossContainer& edge_p_cross, const PInternalContainer& gU_p_internal,
                         size_t gU_n_members, size_t gV_n_members, size_t j0) -> edge_probability_t {
  auto to_p_to_target = TRANSFORM_VIEW(get_merged_p_simple_to(edge_p_cross, gU_p_internal, gU_n_members, j0, _1));
  return at_least_1_probability_r(range(gV_n_members) | to_p_to_target);
}

// The probability of spread success from gu.members[j0] to any member of gv
template <is_coarsened_edge_details EdgeDetailsType, is_coarsened_vertex_details VertexDetailsType>
auto get_merged_p_simple(const EdgeDetailsType& edge, const VertexDetailsType& gu, const VertexDetailsType& gv,
                         size_t j0) -> edge_probability_t {
  return get_merged_p_simple(edge.p_cross, gu.p_internal, gu.n_members(), gv.n_members(), j0);
}

// The probability of spread success from gu.members[j0] to ANY member of gv
// is_seed: Whether the source vertex gu.members[j0] is a seed.
auto get_wim_merged_p_simple(const WIMCoarsenedEdgeDetails& edge, const WIMCoarsenedVertexDetails& gu,
                             const WIMCoarsenedVertexDetails& gv, size_t j0, bool is_seed) -> edge_probability_t {
  auto [p_cross_first, p_internal_first] =
      is_seed ? std::tie(edge.p_seed_cross, gu.p_seed_internal) : std::tie(edge.p_cross, gu.p_internal);
  return get_merged_p_simple_with_first( //
      edge.p_cross, p_cross_first, gu.p_internal, p_internal_first, gu.n_members(), gv.n_members(), j0);
}

// The probability of spread success from gu.members[j0] to ANY member of gv
// i_boosted: Which member in Gv is assumed to be boosted.
auto get_wbim_merged_p_simple(const WBIMCoarsenedEdgeDetails& edge, const WBIMCoarsenedVertexDetails& gu,
                              const WBIMCoarsenedVertexDetails& gv, size_t j0, size_t i_boosted = -1)
    -> edge_probability_t {
  // We only care about p_cross[*][i] or p_boost_cross[*][i]
  auto to_p_to_target = views::transform([&](size_t i) {
    const auto* p_cross_used = (i == i_boosted) ? &edge.p_boost_cross : &edge.p_cross;
    return get_merged_p_simple_to(*p_cross_used, gu.p_internal, gu.n_members(), j0, i);
  });
  return at_least_1_probability_r(gv.member_indices() | to_p_to_target); // To ANY member of Gv
}

// F[2^N][2^N]
using PreciseStateContainer = Array2D<edge_probability_t, (1 << MAX_N_MEMBERS), (1 << MAX_N_MEMBERS)>;

auto make_initial_precise_state_container() -> PreciseStateContainer {
  auto res = PreciseStateContainer{};
  ranges::fill(res | views::join, P_NOT_ASSIGNED);
  return res;
}

struct GetMergedPContext {
  const PCrossContainer* edge_p_cross;
  const PInternalContainer* gU_p_internal;
  size_t gU_n_members;
  size_t gV_n_members;
};

// T = Subset of gu.members that have received message in all the previous turns
// S = Subset of gu.members that receives message in the last turn.
// Both T, S are encoded as bitset. S is always a non-empty subset of T.
auto get_merged_p_precise_impl(const GetMergedPContext* ctx, //
                               PreciseStateContainer& F, size_t T, size_t S) -> edge_probability_t {
  auto N = ctx->gU_n_members;
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
    for (auto i : range(ctx->gV_n_members)) {
      paths.push_back((*ctx->edge_p_cross)[j][i]);
    }
  }
  auto p_success_in_this_turn = at_least_1_probability_r(paths);
  constexpr auto msg_pattern = "F[{1:#0{0}b}][{2:#0{0}b}]: Probability of this turn = {3:.4f}";
  // 2 : extra width of the prefix "0b"
  MYLOG_FMT_TRACE(msg_pattern, N + 2, T, S, p_success_in_this_turn);

  auto sum = 0.0_ep;
  auto to_p_internal = views::transform([&](auto j_pair) {
    auto [j0, j1] = j_pair;
    return (*ctx->gU_p_internal)[j0][j1];
  });
  // Traverses all the non-empty subset of U - T
  for (auto S_next = U ^ T; S_next != 0; S_next = (S_next - 1) & (U ^ T)) {
    auto others = U ^ T ^ S_next;
    // For each dest in S_next, at least 1 path in S -> {dest} must pass
    auto p_S_to_S_next = 1.0_ep;
    for (auto dest : indices_in_bitset(N, S_next)) {
      p_S_to_S_next *= at_least_1_probability_r( //
          indices_in_bitset(N, S) | TRANSFORM_VIEW((*ctx->gU_p_internal)[_1][dest]));
    }
    // Each direct path in S -> others must fail
    auto p_S_to_others = at_least_1_probability_r(
        views::cartesian_product(indices_in_bitset(N, S), indices_in_bitset(N, others)) | to_p_internal);
    // Message propagates to S_next in the next turn
    auto p_next = get_merged_p_precise_impl(ctx, F, T | S_next, S_next);
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

auto get_merged_p_precise(const GetMergedPContext& ctx, PreciseStateContainer* F, size_t T, size_t S)
    -> edge_probability_t {
  if (F == nullptr) {
    auto F_temp = make_initial_precise_state_container();
    return get_merged_p_precise_impl(&ctx, F_temp, T, S);
  }
  return get_merged_p_precise_impl(&ctx, *F, T, S);
}

auto get_merged_p_precise(const GetMergedPContext& ctx, PreciseStateContainer* F, size_t S) -> edge_probability_t {
  return get_merged_p_precise(ctx, F, S, S);
}

auto get_merged_p(const PCrossContainer& edge_p_cross, const PInternalContainer& gU_p_internal,
                  const NormalizedHeuristics& heuristics_in, size_t gU_n_members, size_t gV_n_members,
                  const CoarseningParams& params) -> edge_probability_t {
  auto N = gU_n_members;
  const auto& [p_in, p_in_all_zero] = heuristics_in;

  // If no in-edge to u other than those from Gv, then there's no non-seed path from Gu to Gv
  if (p_in_all_zero) {
    return 0.0_ep;
  }
  // Special case: Gu has only 1 member {u}
  if (N == 1) {
    return at_least_1_probability_r(edge_p_cross[0] | views::take(gV_n_members));
  }
  BOOST_ASSERT_MSG(N == 2 || N == 3, "Group size other than 1, 2, 3 is not supported.");

  auto precise_ctx = GetMergedPContext{
      .edge_p_cross = &edge_p_cross,
      .gU_p_internal = &gU_p_internal,
      .gU_n_members = gU_n_members,
      .gV_n_members = gV_n_members,
  };
  // SEPARATE_SIMPLE, MERGED_SIMPLE, SEPARATE_PRECISE
  if (params.edge_weight_rule != EdgeWeightRule::MERGED_PRECISE) {
    // p_paths[j0] = Probability from Gu to Gv.members[j0]
    auto p_paths = [&]() -> StaticVector<edge_probability_t, MAX_N_MEMBERS> {
      if (params.edge_weight_rule == EdgeWeightRule::SEPARATE_PRECISE) {
        auto F = make_initial_precise_state_container();
        auto view = range(N) | TRANSFORM_VIEW(get_merged_p_precise(precise_ctx, &F, make_bitset_from_indices({_1})));
        return {view.begin(), view.end()};
      } else {
        auto to_merged_p = LAMBDA_1(get_merged_p_simple(edge_p_cross, gU_p_internal, gU_n_members, gV_n_members, _1));
        auto view = range(N) | views::transform(to_merged_p);
        return {view.begin(), view.end()};
      }
    }();

    auto to_combined = TRANSFORM_VIEW(p_in[_1] * p_paths[_1]);
    if (params.edge_weight_rule == EdgeWeightRule::MERGED_SIMPLE) {
      // Conditional probability with virtual source: P(VS -> Gu -> Gv) / P(VS -> Gu)
      return at_least_1_probability_r(range(N) | to_combined) / at_least_1_probability_r(p_in);
    } else {
      // Weighted sum with p_in[i] as weight of gu.members[i]
      return accumulate_sum(range(N) | to_combined);
    }
  } // if (params.edge_weight_rule != EdgeWeightRule::MERGED_PRECISE)

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
    auto p_pass = get_merged_p_precise(precise_ctx, &F, S);
    constexpr auto msg_pattern = "S = {1:#0{0}b}:"
                                 "\n\tp(VS ---> S)      = {2:.4f}"
                                 "\n\tp(VS -/-> Gu - S) = {3:.4f}"
                                 "\n\tp(S  ---> Gv)     = {4:.4f}";
    MYLOG_FMT_TRACE(msg_pattern, N + 2, S, p_enter_1, p_enter_0, p_pass); // 2 : extra width of the prefix "0b"
    sum += p_enter_1 * p_enter_0 * p_pass;
  }
  BOOST_ASSERT_MSG(sum >= 0.0 && sum <= 1.0, "Implementation error: probability sum is out of range [0, 1].");
  // Conditional probability: P(VS -> gu.members -> gv.members) / P(VS -> gu.members)
  return sum / at_least_1_probability_r(p_in);
}

auto merge_coarsened_wim_edge_p(WIMCoarsenedEdgeDetails& dest, const WIMCoarsenedEdgeDetails* inv_dest,
                                const WIMCoarsenedVertexDetails& gu, const WIMCoarsenedVertexDetails& gv,
                                const CoarseningParams& params) -> void {
  auto heuristics = get_normalized_heuristic_p_in( //
      dest, inv_dest, gu, gv, params.in_out_heuristic_rule, NormalizedHeuristics::KEEPS_WHEN_ALL_ZERO);
  dest.merged.p = get_merged_p(dest.p_cross, gu.p_internal, heuristics, gu.n_members(), gv.n_members(), params);
}

auto merge_coarsened_wbim_edge_p(WBIMCoarsenedEdgeDetails& dest, const WBIMCoarsenedEdgeDetails* inv_dest,
                                 const WBIMCoarsenedVertexDetails& gu, const WBIMCoarsenedVertexDetails& gv,
                                 const CoarseningParams& params) -> void {
  auto heuristics = [&]() {
    if (gu.is_seed || rule_prefix_is_seed(params.in_out_heuristic_rule)) {
      return get_normalized_heuristic_p_in( //
          gu, params.in_out_heuristic_rule, NormalizedHeuristics::TO_UNIFORM_WHEN_ALL_ZERO);
    }
    return get_normalized_heuristic_p_in( //
        dest, inv_dest, gu, gv, params.in_out_heuristic_rule, NormalizedHeuristics::KEEPS_WHEN_ALL_ZERO);
  }();
  BOOST_ASSERT_MSG(!gu.is_seed || (heuristics.p.size() == 1 && heuristics.p[0] == 1.0_ep),
                   "Seed group with more than 1 member is disallowed.");
  dest.merged.p = get_merged_p(dest.p_cross, gu.p_internal, heuristics, gu.n_members(), gv.n_members(), params);
}

auto merge_coarsened_wim_edge_p_seed(WIMCoarsenedEdgeDetails& dest, const WIMCoarsenedVertexDetails& gu,
                                     const WIMCoarsenedVertexDetails& gv, const CoarseningParams& params) -> void {
  auto N = gu.n_members();
  // Special case: Gu has only 1 member {u}
  if (N == 1) {
    // Note: the seed path is always enabled, whether or not there's in-edge to u.
    dest.merged.p_seed = at_least_1_probability_r(dest.p_seed_cross[0] | views::take(gv.n_members()));
    return;
  }
  BOOST_ASSERT_MSG(N == 2 || N == 3, "Group size other than 1, 2, 3 is not supported.");
  dest.merged.p_seed = [&]() {
    if (params.seed_edge_weight_rule == SeedEdgeWeightRule::BEST_SEED_INDEX) {
      auto res = get_wim_merged_p_simple(dest, gu, gv, gu.best_seed_index, true);
      MYLOG_FMT_TRACE("p_seed from group {} to {} = {:.4f}, starting from the locally best seed candidate #{}", //
                      gu.members, gv.members, res, gu.members[gu.best_seed_index]);
      return res;
    }
    auto estimated_gain_as_seed = StaticVector<edge_probability_t, MAX_N_MEMBERS>(N);
    for (auto j0 : range(N)) {
      estimated_gain_as_seed[j0] = get_wim_merged_p_simple(dest, gu, gv, j0, true);
    }
    MYLOG_FMT_TRACE("Estimated gain as seed of group {} = {}", gu.members, estimated_gain_as_seed);
    if (params.seed_edge_weight_rule == SeedEdgeWeightRule::AVERAGE) {
      return accumulate_sum(estimated_gain_as_seed) / N;
    }
    BOOST_ASSERT(params.seed_edge_weight_rule == SeedEdgeWeightRule::MAX);
    return ranges::max(estimated_gain_as_seed);
  }();
}

auto merge_coarsened_wbim_edge_p_boost(WBIMCoarsenedEdgeDetails& dest, const WBIMCoarsenedEdgeDetails* inv_dest,
                                       const WBIMCoarsenedVertexDetails& gu, const WBIMCoarsenedVertexDetails& gv,
                                       const CoarseningParams& params) -> void {
  auto N = gu.n_members();
  BOOST_ASSERT_MSG(N >= 1 && N <= 3, "Group size other than 1, 2, 3 is not supported.");
  // Assumes the Gv.members[i] is boosted
  auto heuristics = [&]() {
    if (gu.is_seed || rule_prefix_is_seed(params.in_out_heuristic_rule)) {
      return get_normalized_heuristic_p_in( //
          gu, params.in_out_heuristic_rule, NormalizedHeuristics::TO_UNIFORM_WHEN_ALL_ZERO);
    }
    return get_normalized_heuristic_p_in( //
        dest, inv_dest, gu, gv, params.in_out_heuristic_rule, NormalizedHeuristics::KEEPS_WHEN_ALL_ZERO);
  }();
  BOOST_ASSERT_MSG(!gu.is_seed || (heuristics.p.size() == 1 && heuristics.p[0] == 1.0_ep),
                   "Seed group with more than 1 member is disallowed.");

  auto get_p_boost = [&](size_t i) {
    auto boosted_p_cross = dest.p_cross;
    ranges::for_each(range(N), LAMBDA_1(boosted_p_cross[_1][i] = dest.p_boost_cross[_1][i]));
    MYLOG_FMT_TRACE("p_cross from group {} to {} when target vertex {} is assumed to be boosted:\n\t{:::.4f}", //
                    gu.members, gv.members, gv.members[i], boosted_p_cross);
    return get_merged_p(boosted_p_cross, gu.p_internal, heuristics, gu.n_members(), gv.n_members(), params);
  };
  dest.merged.p_boost = [&]() {
    if (params.boosted_edge_weight_rule == BoostedEdgeWeightRule::BEST_BOOSTED_INDEX) {
      return get_p_boost(gv.best_boosted_index);
    }
    // Attempts each boosted target in Gv
    auto p_as_boosted = StaticVector<edge_probability_t, MAX_N_MEMBERS>(gv.n_members());
    ranges::for_each(gv.member_indices(), LAMBDA_1(p_as_boosted[_1] = get_p_boost(_1)));
    MYLOG_FMT_TRACE("Estimated p_boost from group {} to {}: {::.4f}", gu.members, gv.members, p_as_boosted);

    if (params.boosted_edge_weight_rule == BoostedEdgeWeightRule::AVERAGE) {
      return accumulate_sum(p_as_boosted) / gv.n_members();
    } else {
      BOOST_ASSERT(params.boosted_edge_weight_rule == BoostedEdgeWeightRule::MAX);
      return ranges::max(p_as_boosted);
    }
  }();
}

auto merge_coarsened_wim_graph_edge(WIMCoarsenedEdgeDetails& dest, const WIMCoarsenedEdgeDetails* inv_dest,
                                    const WIMCoarsenedVertexDetails& gu, const WIMCoarsenedVertexDetails& gv,
                                    const CoarseningParams& params) -> void {
  // Part (1) merged.p, see above
  merge_coarsened_wim_edge_p(dest, inv_dest, gu, gv, params);
  // Part (2): merged.p_seed, see above
  merge_coarsened_wim_edge_p_seed(dest, gu, gv, params);
  // Part (3): Ensures 0.0 <= p <= p_seed <= 1.0
  dest.merged.p = std::clamp(dest.merged.p, 0.0_ep, 1.0_ep);
  dest.merged.p_seed = std::clamp(dest.merged.p_seed, dest.merged.p, 1.0_ep);
}

auto merge_coarsened_wbim_graph_edge(WBIMCoarsenedEdgeDetails& dest, const WBIMCoarsenedEdgeDetails* inv_dest,
                                     const WBIMCoarsenedVertexDetails& gu, const WBIMCoarsenedVertexDetails& gv,
                                     const CoarseningParams& params) -> void {
  // Part (1) merged.p, see above
  merge_coarsened_wbim_edge_p(dest, inv_dest, gu, gv, params);
  // Part (2): merged.p_seed, see above
  merge_coarsened_wbim_edge_p_boost(dest, inv_dest, gu, gv, params);
  // Part (3): Ensures 0.0 <= p <= p_boost <= 1.0
  dest.merged.p = std::clamp(dest.merged.p, 0.0_ep, 1.0_ep);
  dest.merged.p_boost = std::clamp(dest.merged.p_boost, dest.merged.p, 1.0_ep);
}

// Dispatching
template <is_coarsened_edge_details EdgeDetailsType, is_coarsened_vertex_details VertexDetailsType>
auto merge_coarsened_graph_edge(EdgeDetailsType& dest, const EdgeDetailsType* inv_dest, const VertexDetailsType& gu,
                                const VertexDetailsType& gv, const CoarseningParams& params) -> void {
  if constexpr (std::is_same_v<VertexDetailsType, WIMCoarsenedVertexDetails>) {
    return merge_coarsened_wim_graph_edge(dest, inv_dest, gu, gv, params);
  } else if constexpr (std::is_same_v<VertexDetailsType, WBIMCoarsenedVertexDetails>) {
    return merge_coarsened_wbim_graph_edge(dest, inv_dest, gu, gv, params);
  } else {
    static_assert(rfl::always_false_v<VertexDetailsType>, "Invalid vertex details type.");
  }
}

// If is_inverse, then the edge source -> target is mapped to Gv -> Gu
// Otherwise, the edge source -> target is mapped to Gu -> Gv
template <is_edge_property E>
struct EdgeCoarseningItemInfo {
  bool is_inverse;
  vertex_id_t source;
  vertex_id_t target;
  E w;

  auto dump() const {
    constexpr auto msg_pattern = "{{.is_inverse = {}, .source = {}, .target = {}, .w = {}}}";
    return fmt::format(msg_pattern, is_inverse, source, target, w);
  }
};

template <is_edge_property E, same_as_either<std::nullptr_t, VertexSet> TargetSkipSet = std::nullptr_t>
auto make_sorted_edge_items(AdjacencyList<E>& graph, vertex_id_t n_coarsened, std::span<const vertex_id_t> group_id,
                            const TargetSkipSet& target_skip_set = nullptr) {
  using ItemInfo = EdgeCoarseningItemInfo<E>;
  using ItemType = VertexPairWithValue<const ItemInfo*>;

  auto [n, m] = graph_n_m(graph);
  BOOST_ASSERT_MSG(group_id.size() == n, "Size mismatch between |V| and group_id list.");
  if constexpr (std::is_same_v<TargetSkipSet, VertexSet>) {
    BOOST_ASSERT_MSG(target_skip_set.n_vertices_in_whole_graph() == n,
                     "Size mismatch between |V| and target_skip_set.");
  }

  // Note: RESERVATION is REQUIRED to prevent pointer invalidation
  auto item_info_values = make_reserved_vector<ItemInfo>(m);
  auto items = make_reserved_vector<ItemType>(m);

  for (auto [u, v, w] : graph::make_edge_range<0>(graph)) {
    auto gu = group_id[u], gv = group_id[v];
    if (gu == gv) {
      continue; // Skips the in-group edges
    }
    if constexpr (std::is_same_v<TargetSkipSet, VertexSet>) {
      if (target_skip_set.contains(v)) {
        continue; // Skips the edges whose target is specified as skipped
      }
    }
    if (gu < gv) {
      auto info = ItemInfo{.is_inverse = false, .source = u, .target = v, .w = w};
      auto pos = &item_info_values.emplace_back(info);
      items.push_back(ItemType{.u = gu, .v = gv, .value = pos});
    } else {
      auto info = ItemInfo{.is_inverse = true, .source = u, .target = v, .w = w};
      auto pos = &item_info_values.emplace_back(info);
      items.push_back(ItemType{.u = gv, .v = gu, .value = pos});
    }
  }
  radix_sort_vertex_pairs(items, n_coarsened);
  MYLOG_TRACE([&]() {
    auto to_item_str = LAMBDA_1(fmt::format("\n\t{{.u = {}, .v = {}, .value -> {}}}", _1.u, _1.v, _1.value->dump()));
    auto contents = items | views::transform(to_item_str) | views::join | views::common;
    return "Items after radix sorting: ["s + std::string(contents.begin(), contents.end()) + "\n]";
  }());
  // The internal data pointer is unchanged after moving
  return std::tuple{std::move(item_info_values), std::move(items)};
}

template <is_edge_property E, class CoarseningDetailsType,
          same_as_either<std::nullptr_t, VertexSet> TargetSkipSet = std::nullptr_t>
auto get_coarsened_graph_edges_impl(AdjacencyList<E>& graph, const CoarseningDetailsType& details,
                                    const CoarseningParams& params, const TargetSkipSet& target_skip_set = nullptr)
    -> DirectedEdgeList<E> {
  using EdgeDetailsType = typename CoarseningDataTraits<E>::EdgeDetails;
  using ItemInfo = EdgeCoarseningItemInfo<E>;
  using ItemType = VertexPairWithValue<const ItemInfo*>;
  auto [item_info_values, items] = make_sorted_edge_items( //
      graph, details.n_coarsened, details.group_id, target_skip_set);

  auto set_p_cross = [&](EdgeDetailsType& dest, const E& w, size_t j_source, size_t j_target) {
    if constexpr (requires { dest.p_cross; }) {
      dest.p_cross[j_source][j_target] = w.p;
    }
    if constexpr (requires { dest.p_seed_cross; }) {
      dest.p_seed_cross[j_source][j_target] = w.p_seed;
    }
    if constexpr (requires { dest.p_boost_cross; }) {
      dest.p_boost_cross[j_source][j_target] = w.p_boost;
    }
  };

  auto res = DirectedEdgeList<E>{details.n_coarsened}; // Keeps coarsened |V| unchanged
  res.open_for_push_back();
  // Each chunk represents all the edges of Gu -> Gv (is_inverse == false) and Gv -> Gu (is_inverse == true)
  for (auto chunk : items | CHUNK_BY_VIEW(_1.u == _2.u && _1.v == _2.v)) {
    auto [gu, gv] = chunk[0].vertices();
    auto [Mu, Mv] = details.group_id_to_size(gu, gv);
    auto E_uv = EdgeDetailsType(Mu, Mv); // Gu -> Gv
    auto E_vu = EdgeDetailsType(Mv, Mu); // Gv -> Gu

    auto has_uv = false;
    auto has_vu = false;
    for (auto info : chunk | views::transform(&ItemType::value)) {
      auto [j_source, j_target] = details.to_index_in_group(info->source, info->target);
      if (info->is_inverse) {
        set_p_cross(E_vu, info->w, j_source, j_target); // The edge source -> target is mapped to Gv -> Gu
        has_vu = true;
      } else {
        set_p_cross(E_uv, info->w, j_source, j_target); // The edge source -> target is mapped to Gu -> Gv
        has_uv = true;
      }
    }
    if (has_uv) {
      merge_coarsened_graph_edge(E_uv, &E_vu, details.groups[gu], details.groups[gv], params);
      res.push_back(gu, gv, E_uv.merged);
      MYLOG_FMT_TRACE("Details of coarsened edge {} -> {}: {:4}", gu, gv, E_uv);
    }
    if (has_vu) {
      merge_coarsened_graph_edge(E_vu, &E_uv, details.groups[gv], details.groups[gu], params);
      res.push_back(gv, gu, E_vu.merged);
      MYLOG_FMT_TRACE("Details of coarsened edge {} -> {}: {:4}", gv, gu, E_vu);
    }
  }
  res.close_for_push_back();
  return res;
}

auto get_coarsened_wim_graph_edges(AdjacencyList<WIMEdge>& graph, const WIMCoarseningDetails& details,
                                   const CoarseningParams& params) -> DirectedEdgeList<WIMEdge> {
  return get_coarsened_graph_edges_impl(graph, details, params);
}

auto get_coarsened_wbim_graph_edges(AdjacencyList<WBIMEdge>& graph, const WBIMCoarseningDetails& details,
                                    const VertexSet& seeds, const CoarseningParams& params)
    -> DirectedEdgeList<WBIMEdge> {
  // Edges whose target is a seed vertex are skipped since they have no influence.
  return get_coarsened_graph_edges_impl(graph, details, params, seeds);
}

template <is_coarsened_vertex_details VertexDetailsType>
auto get_coarsened_vertex_weights(const CoarseningDetails<VertexDetailsType>& details,
                                  std::span<const vertex_weight_t> vertex_weights, const CoarseningParams& params,
                                  NormalizedHeuristics::PolicyWhenAllZero policy) -> std::vector<vertex_weight_t> {
  auto coarsened_vertex_weights = std::vector<vertex_weight_t>(details.n_coarsened, 0.0);
  // Traverses each group. vc = Coarsened vertex index in range [0, Nc)
  for (auto [vc, gc] : details.groups | views::enumerate) {
    auto& dest = coarsened_vertex_weights[vc];
    auto [p_in, p_in_all_zero] = get_normalized_heuristic_p_in(gc, params.in_out_heuristic_rule, policy);

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

struct CoarsenWIMGraphImplResult {
  DirectedEdgeList<WIMEdge> coarsened_edge_list;
  std::vector<vertex_weight_t> coarsened_vertex_weights;
  WIMCoarseningDetails details;
};

auto coarsen_wim_graph_impl(AdjacencyList<WIMEdge>& graph, InvAdjacencyList<WIMEdge>& inv_graph,
                            std::span<const vertex_weight_t> vertex_weights, vertex_id_t n_groups,
                            std::span<const vertex_id_t> group_id, const CoarseningParams& params)
    -> rfl::Result<CoarsenWIMGraphImplResult> {
  MYLOG_FMT_DEBUG("Parameters: {:4}", params);
  if (rule_prefix_is_seed(params.in_out_heuristic_rule)) {
    return rfl::Error{"InOutHeuristicRule::SEED_* is disallowed in WIM graph coarsening "
                      "since there's no seed known in prior."};
  }
  auto n = graph::num_vertices(graph);
  auto details = WIMCoarseningDetails(               //
      WIMCoarseningDetails::value_initialize_groups, //
      n, n_groups, std::vector(group_id.begin(), group_id.end()), std::vector<vertex_id_t>(n));
  auto step_timer = StepTimer("WIM coarsening after matching");
  MYLOG_FMT_DEBUG("Starts WIM coarsening: n = {}, n_groups = {}", details.n, details.n_coarsened);

  // Note: For unweighted graph, it's required to create a dummy list like std::vector<vertex_weight_t>(n, 1.0_vw)
  BOOST_ASSERT_MSG(vertex_weights.size() == n, "Size mismatch");

  // Step 1: groups[].members & groups[].vertex_weights & groups[].heuristics_* & index_in_group
  step_timer.do_step("assigning group members", [&]() {
    // Vertex v is placed in group g
    for (auto [v, g] : views::enumerate(group_id)) {
      auto& cur_group = details.groups[g];
      details.index_in_group[v] = static_cast<vertex_id_t>(cur_group.n_members());
      cur_group.members.push_back(v);
      cur_group.vertex_weights.push_back(vertex_weights[v]);
      // Heuristics
      cur_group.heuristics_in.push_back( //
          get_heuristics_in(inv_graph, group_id, vertex_weights, v, params.in_out_heuristic_rule));
      cur_group.heuristics_out.push_back( //
          get_heuristics_out(graph, group_id, vertex_weights, v, params.in_out_heuristic_rule));
      cur_group.heuristics_out_seed.push_back( //
          get_heuristics_out_seed(graph, group_id, vertex_weights, v, params.in_out_heuristic_rule));
    }
  });
  // Step 2: groups[].p_internal
  step_timer.do_step("assigning p_internal of groups", [&]() {
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
  step_timer.do_step("selecting best seed candidate in each group", [&]() {
    for (auto& g : details.groups) {
      g.best_seed_index = select_best_seed_in_group(g).index_in_group;
    }
  });
  // Step 4: Generates coarsened vertex weights
  auto coarsened_vertex_weights = step_timer.do_step("coarsening vertex weights", [&]() {
    return get_coarsened_vertex_weights( //
        details, vertex_weights, params, NormalizedHeuristics::KEEPS_WHEN_ALL_ZERO);
  });
  // Step 5: Finally, builds the result from all the components
  auto coarsened_edge_list = step_timer.do_step("building edge list of the coarsened graph", [&]() {
    return get_coarsened_wim_graph_edges(graph, details, params);
  });

  return CoarsenWIMGraphImplResult{
      .coarsened_edge_list = std::move(coarsened_edge_list),
      .coarsened_vertex_weights = std::move(coarsened_vertex_weights),
      .details = std::move(details),
  };
}

template <class ResultType>
auto coarsen_wim_graph(AdjacencyList<WIMEdge>& graph, InvAdjacencyList<WIMEdge>& inv_graph,
                       std::span<const vertex_weight_t> vertex_weights, vertex_id_t n_groups,
                       std::span<const vertex_id_t> group_id, const CoarseningParams& params)
    -> rfl::Result<ResultType> {
  return coarsen_wim_graph_impl(graph, inv_graph, vertex_weights, n_groups, group_id, params)
      .transform([&](CoarsenWIMGraphImplResult res) {
        auto res_details = [&]() {
          if constexpr (std::is_same_v<ResultType, WIMCoarsenGraphDetailedResult>) {
            return std::move(res.details); // Returns the details directly
          } else if constexpr (std::is_same_v<ResultType, WIMCoarsenGraphBriefResult>) {
            return res.details.to_brief(); // Transforms detailed to brief information
          } else {
            static_assert(rfl::always_false_v<ResultType>, "Invalid result type.");
          }
        }();
        return ResultType(res.coarsened_edge_list, std::move(res.coarsened_vertex_weights), std::move(res_details));
      });
}

struct CoarsenWBIMGraphImplResult {
  DirectedEdgeList<WBIMEdge> coarsened_edge_list;
  std::vector<vertex_weight_t> coarsened_vertex_weights;
  VertexSet coarsened_seeds;
  WBIMCoarseningDetails details;
};

auto coarsen_wbim_graph_impl(AdjacencyList<WBIMEdge>& graph, InvAdjacencyList<WBIMEdge>& inv_graph,
                             std::span<const vertex_weight_t> vertex_weights, const VertexSet& seeds,
                             vertex_id_t n_groups, std::span<const vertex_id_t> group_id,
                             const CoarseningParams& params) -> rfl::Result<CoarsenWBIMGraphImplResult> {
  auto n = graph::num_vertices(graph);
  auto details = WBIMCoarseningDetails(               //
      WBIMCoarseningDetails::value_initialize_groups, //
      n, n_groups, std::vector(group_id.begin(), group_id.end()), std::vector<vertex_id_t>(n));
  auto step_timer = StepTimer("WBIM coarsening after matching");
  MYLOG_FMT_DEBUG("Starts WBIM coarsening: n = {}, n_groups = {}, parameters: {:4}", //
                  details.n, details.n_coarsened, params);

  // Note: For unweighted graph, it's required to create a dummy list like std::vector<vertex_weight_t>(n, 1.0_vw)
  BOOST_ASSERT_MSG(vertex_weights.size() == n, "Size mismatch");

  // Step 1: groups[].members & groups[].vertex_weights & groups[].heuristics_* & index_in_group
  auto step_1_res = step_timer.do_step("assigning group members", [&]() -> ResultVoid {
    // For SEED rule: Calculates heuristic by seed detection
    auto activation = [&]() -> rfl::Result<WBIMActivationProbability> {
      if (rule_prefix_is_seed(params.in_out_heuristic_rule)) {
        return wbim_activation_probability_from_seeds(inv_graph, seeds, params.max_distance_from_seed);
      }
      return WBIMActivationProbability{}; // Dummy if unused
    }();
    return activation.transform([&](const WBIMActivationProbability& activation) {
      // Vertex v is placed in group g
      for (auto [v, g] : views::enumerate(group_id)) {
        auto& cur_group = details.groups[g];
        details.index_in_group[v] = static_cast<vertex_id_t>(cur_group.n_members());
        cur_group.members.push_back(v);
        cur_group.vertex_weights.push_back(vertex_weights[v]);
        // Heuristics
        if (rule_prefix_is_seed(params.in_out_heuristic_rule)) {
          cur_group.heuristics_in.push_back(activation.p_in[v]);
          cur_group.heuristics_in_boost.push_back(activation.p_in_boosted[v]);

          auto rule_out = [&params]() {
            if (params.in_out_heuristic_rule == InOutHeuristicRule::SEED_P) {
              return InOutHeuristicRule::P;
            }
            BOOST_ASSERT(params.in_out_heuristic_rule == InOutHeuristicRule::SEED_W);
            return InOutHeuristicRule::W;
          }();
          // Seeds as out-neighbors are skipped for better heuristics accuracy
          cur_group.heuristics_out.push_back( //
              get_heuristics_out(graph, group_id, vertex_weights, seeds, v, rule_out));
        } else {
          cur_group.heuristics_in.push_back( //
              get_heuristics_in(inv_graph, group_id, vertex_weights, v, params.in_out_heuristic_rule));
          cur_group.heuristics_in_boost.push_back( //
              get_heuristics_in_boost(inv_graph, group_id, vertex_weights, v, params.in_out_heuristic_rule));
          // Seeds as out-neighbors are skipped for better heuristics accuracy
          cur_group.heuristics_out.push_back( //
              get_heuristics_out(graph, group_id, vertex_weights, seeds, v, params.in_out_heuristic_rule));
        }
      }
      return RESULT_VOID_SUCCESS;
    });
  });
  RFL_RETURN_ON_ERROR(step_1_res);

  // Step 2: groups[].p_internal
  step_timer.do_step("assigning p_internal of groups", [&]() {
    for (auto [u, v, w] : graph::make_edge_range<0>(graph)) {
      auto g = group_id[u];
      if (g != group_id[v]) {
        continue;
      }
      auto [ju, jv] = details.to_index_in_group(u, v);
      details.groups[g].p_internal[ju][jv] = w.p;
      details.groups[g].p_boost_internal[ju][jv] = w.p_boost;
    }
  });
  // Step 3: coarsening seeds
  auto coarsened_seeds = step_timer.do_step("mapping seeds to coarsened index", [&]() {
    auto coarsened_seeds = make_reserved_vector<vertex_id_t>(n_groups);
    for (auto [group_index, g] : details.groups | views::enumerate) {
      BOOST_ASSERT_MSG(g.n_members() > 0, "Why an empty group here?");
      if (seeds.contains(g.members[0])) {
        BOOST_ASSERT_MSG(g.n_members() == 1, "Seeds must be grouped alone.");
        coarsened_seeds.push_back(group_index);
        g.is_seed = true;
      } else {
        BOOST_ASSERT_MSG(ranges::none_of(g.members, LAMBDA_1(seeds.contains(_1))), "Seeds must be grouped alone.");
        g.is_seed = false;
      }
    }
    MYLOG_FMT_DEBUG("Coarsened seeds: {}", coarsened_seeds);
    return VertexSet(n_groups, std::move(coarsened_seeds)); // Coarsened |V| = n_groups
  });
  // Step 4: Calculates the best boosted vertex locally in each group
  step_timer.do_step("selecting best boosted candidate in each group", [&]() {
    for (auto& g : details.groups) {
      g.best_boosted_index = select_best_boosted_in_group(g, params.boosted_selection_rule).index_in_group;
    }
  });
  // Step 5: Generates coarsened vertex weights
  auto coarsened_vertex_weights = step_timer.do_step("coarsening vertex weights", [&]() {
    auto policy = rule_prefix_is_seed(params.in_out_heuristic_rule) //
                      ? NormalizedHeuristics::TO_UNIFORM_WHEN_ALL_ZERO
                      : NormalizedHeuristics::KEEPS_WHEN_ALL_ZERO;
    return get_coarsened_vertex_weights(details, vertex_weights, params, policy);
  });
  // Step 6: Finally, builds the result from all the components
  auto coarsened_edge_list = step_timer.do_step("building edge list of the coarsened graph", [&]() {
    return get_coarsened_wbim_graph_edges(graph, details, seeds, params); // Skips the edges whose target is a seed
  });

  return CoarsenWBIMGraphImplResult{
      .coarsened_edge_list = std::move(coarsened_edge_list),
      .coarsened_vertex_weights = std::move(coarsened_vertex_weights),
      .coarsened_seeds = std::move(coarsened_seeds),
      .details = std::move(details),
  };
}

template <class ResultType>
auto coarsen_wbim_graph(AdjacencyList<WBIMEdge>& graph, InvAdjacencyList<WBIMEdge>& inv_graph,
                        std::span<const vertex_weight_t> vertex_weights, const VertexSet& seeds, vertex_id_t n_groups,
                        std::span<const vertex_id_t> group_id, const CoarseningParams& params)
    -> rfl::Result<ResultType> {
  return coarsen_wbim_graph_impl(graph, inv_graph, vertex_weights, seeds, n_groups, group_id, params)
      .transform([&](CoarsenWBIMGraphImplResult res) {
        auto res_details = [&]() {
          if constexpr (std::is_same_v<ResultType, WBIMCoarsenGraphDetailedResult>) {
            return std::move(res.details); // Returns the details directly
          } else if constexpr (std::is_same_v<ResultType, WBIMCoarsenGraphBriefResult>) {
            return res.details.to_brief(); // Transforms detailed to brief information
          } else {
            static_assert(rfl::always_false_v<ResultType>, "Invalid result type.");
          }
        }();
        return ResultType(res.coarsened_edge_list, std::move(res.coarsened_vertex_weights),
                          std::move(res.coarsened_seeds), std::move(res_details));
      });
}

template <class ResultType>
auto coarsen_wim_graph_by_match_generic(const AdjacencyList<WIMEdge>& graph, const InvAdjacencyList<WIMEdge>& inv_graph,
                                        std::span<const vertex_weight_t> vertex_weights, vertex_id_t n_groups,
                                        std::span<const vertex_id_t> group_id, const CoarseningParams& params) noexcept
    -> rfl::Result<ResultType> {
  auto n = graph::num_vertices(graph);
  if (vertex_weights.empty()) {
    auto dummy_vertex_weights = std::vector<vertex_weight_t>(n, 1.0_vw);
    return coarsen_wim_graph<ResultType>( // Uses dummy vertex weight range
        as_non_const(graph), as_non_const(inv_graph), dummy_vertex_weights, n_groups, group_id, params);
  } else if (vertex_weights.size() != n) {
    constexpr auto msg_pattern = "Size mismatch between vertex weights (which is {}) and |V| (which is {})";
    return rfl::Error{fmt::format(msg_pattern, vertex_weights.size(), n)};
  }
  return coarsen_wim_graph<ResultType>( //
      as_non_const(graph), as_non_const(inv_graph), vertex_weights, n_groups, group_id, params);
}

template <class ResultType>
auto coarsen_wbim_graph_by_match_generic(const AdjacencyList<WBIMEdge>& graph,
                                         const InvAdjacencyList<WBIMEdge>& inv_graph,
                                         std::span<const vertex_weight_t> vertex_weights, const VertexSet& seeds,
                                         vertex_id_t n_groups, std::span<const vertex_id_t> group_id,
                                         const CoarseningParams& params) noexcept -> rfl::Result<ResultType> {
  auto n = graph::num_vertices(graph);
  if (seeds.n_vertices_in_whole_graph() != n) {
    constexpr auto msg_pattern = "Size mismatch between seeds (which is {}) and |V| (which is {})";
    return rfl::Error{fmt::format(msg_pattern, seeds.n_vertices_in_whole_graph(), n)};
  }
  if (vertex_weights.empty()) {
    auto dummy_vertex_weights = std::vector<vertex_weight_t>(n, 1.0_vw);
    return coarsen_wbim_graph<ResultType>( // Uses dummy vertex weight range
        as_non_const(graph), as_non_const(inv_graph), dummy_vertex_weights, seeds, n_groups, group_id, params);
  } else if (vertex_weights.size() != n) {
    constexpr auto msg_pattern = "Size mismatch between vertex weights (which is {}) and |V| (which is {})";
    return rfl::Error{fmt::format(msg_pattern, vertex_weights.size(), n)};
  }
  return coarsen_wbim_graph<ResultType>( //
      as_non_const(graph), as_non_const(inv_graph), vertex_weights, seeds, n_groups, group_id, params);
}

template <class ResultType>
auto coarsen_wim_graph_framework(const AdjacencyList<WIMEdge>& graph, const InvAdjacencyList<WIMEdge>& inv_graph,
                                 std::span<const vertex_weight_t> vertex_weights, const CoarseningParams& params)
    -> rfl::Result<ResultType> {
  auto [n, m] = graph_n_m(graph);
  auto bidir_graph = merge_wim_edge_to_undirected(graph, params);
  auto [n_groups, group_id] = mongoose_match(bidir_graph, params);
  return coarsen_wim_graph_by_match_generic<ResultType>(graph, inv_graph, vertex_weights, n_groups, group_id, params);
}

template <class ResultType>
auto coarsen_wbim_graph_framework(const AdjacencyList<WBIMEdge>& graph, const InvAdjacencyList<WBIMEdge>& inv_graph,
                                  std::span<const vertex_weight_t> vertex_weights, const VertexSet& seeds,
                                  const CoarseningParams& params) -> rfl::Result<ResultType> {
  auto [n, m] = graph_n_m(graph);
  auto bidir_graph = merge_wbim_edge_to_undirected(graph, params);
  auto [n_groups, group_id] = mongoose_match(bidir_graph, params);
  return coarsen_wbim_graph_by_match_generic<ResultType>( //
      graph, inv_graph, vertex_weights, seeds, n_groups, group_id, params);
}

template <class CoarsenResultType>
auto expand_wim_seed_vertices_impl(const AdjacencyList<WIMEdge>& graph, std::span<const vertex_weight_t> vertex_weights,
                                   const CoarsenResultType& coarsening_result,
                                   std::span<const vertex_id_t> coarsened_seeds, const ExpandingParams& params) noexcept
    -> rfl::Result<ExpandSeedResult> {
  auto n = coarsening_result.details.n;            // # of vertices BEFORE coarsening
  auto nc = coarsening_result.details.n_coarsened; // # of vertices AFTER coarsening
  BOOST_ASSERT_MSG(ranges::size(coarsened_seeds) <= nc, "Expects no more than nc vertices as seeds");
  BOOST_ASSERT_MSG(ranges::max(coarsened_seeds) < nc, "Indices of seeds are out of range [0, nc).");

  auto fallback_vertex_weights = std::vector<vertex_weight_t>{};
  if (vertex_weights.empty()) {
    // If no vertex weight is provided, a list of unit weights is provided
    fallback_vertex_weights.assign(n, 1.0_vw);
    vertex_weights = fallback_vertex_weights;
  } else if (vertex_weights.size() != n) {
    // Raises an error on size mismatch
    constexpr auto msg_pattern = "Size mismatch between vertex weights (which is {}) and |V| (which is {})";
    return rfl::Error{fmt::format(msg_pattern, vertex_weights.size(), n)};
  }

  auto expanded_seeds = make_reserved_vector<vertex_id_t>(coarsened_seeds.size());
  auto make_expanded_seeds_as_s_local = [&]() {
    for (auto sc : coarsened_seeds) {
      auto& group = coarsening_result.details.groups[sc];
      expanded_seeds.push_back(group.members[group.best_seed_index]);
    }
  };

  // RANDOM: Random choice (used for contrast experiment)
  if (params.vertex_expanding_rule == VertexExpandingRule::RANDOM) {
    for (auto sc : coarsened_seeds) {
      auto& group = coarsening_result.details.groups[sc];
      auto rand_index = rand_engine() % group.n_members();
      expanded_seeds.push_back(group.members[rand_index]);
    }
    goto expand_wim_seed_vertices_impl_return;
  }

  // LOCAL: Only considers information inside current group
  if (params.vertex_expanding_rule == VertexExpandingRule::LOCAL) {
    make_expanded_seeds_as_s_local();
    goto expand_wim_seed_vertices_impl_return;
  }

  // SIMULATIVE: Picks one-by-one separately by simulation result of single expanded seed vertex
  if (params.vertex_expanding_rule == VertexExpandingRule::SIMULATIVE) {
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
    goto expand_wim_seed_vertices_impl_return;
  }

  // ITERATIVE: Each iteration attempts to change a vertex
  if (params.vertex_expanding_rule == VertexExpandingRule::ITERATIVE) {
    make_expanded_seeds_as_s_local(); // Result of LOCAL as initial solution
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
      goto expand_wim_seed_vertices_impl_return;
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
      MYLOG_FMT_DEBUG("Iteration #{}: Attempting to change {} to {} (in group {} coarsened to vertex #{})", //
                      iteration_index, prev, next, group.members, coarsened_seeds[si]);
      expanded_seeds[si] = next;
      auto next_sim_res = simulate(graph, vertex_weights, expanded_seeds, nullptr, **params.simulation_try_count);
      RFL_RETURN_ON_ERROR(next_sim_res);
      MYLOG_FMT_DEBUG("Iteration #{}: Candidate {} vs. {}: Simulation result: {:.4f} vs. {:.4f}", //
                      iteration_index, next, prev, *next_sim_res, *sim_res);
      if (*next_sim_res > *sim_res) {
        sim_res = next_sim_res; // Changes prev to next
      } else {
        expanded_seeds[si] = prev; // Otherwise, restores current expanded seed
      }
      MYLOG_FMT_DEBUG("Expanded seed set after iteration #{}: {}", iteration_index, expanded_seeds);
    }
    goto expand_wim_seed_vertices_impl_return;
  }
  BOOST_ASSERT_MSG(false, "Unreachable branch due to invalid expanding rule.");
expand_wim_seed_vertices_impl_return:
  return ExpandSeedResult{.expanded_seeds = std::move(expanded_seeds)};
}

template <class CoarsenResultType>
auto expand_wbim_boosted_vertices_impl(const CoarsenResultType& coarsening_result,
                                       std::span<const vertex_id_t> coarsened_boosted,
                                       const ExpandingParams& params) noexcept -> rfl::Result<ExpandBoostedResult> {
  auto n = coarsening_result.details.n;            // # of vertices BEFORE coarsening
  auto nc = coarsening_result.details.n_coarsened; // # of vertices AFTER coarsening
  BOOST_ASSERT_MSG(ranges::size(coarsened_boosted) <= nc, "Expects no more than nc vertices as boosted vertices");
  BOOST_ASSERT_MSG(ranges::max(coarsened_boosted) < nc, "Indices of boosted vertices are out of range [0, nc).");

  auto expanded_boosted = make_reserved_vector<vertex_id_t>(coarsened_boosted.size());
  if (params.vertex_expanding_rule == VertexExpandingRule::RANDOM) {
    for (auto sc : coarsened_boosted) {
      auto& group = coarsening_result.details.groups[sc];
      auto rand_index = rand_engine() % group.n_members();
      expanded_boosted.push_back(group.members[rand_index]);
    }
  } else if (params.vertex_expanding_rule == VertexExpandingRule::LOCAL) {
    for (auto sc : coarsened_boosted) {
      auto& group = coarsening_result.details.groups[sc];
      expanded_boosted.push_back(group.members[group.best_boosted_index]);
    }
  } else {
    constexpr auto msg_pattern = "Unsupported WBIM expanding rule: {}";
    return rfl::Error{fmt::format(msg_pattern, magic_enum::enum_name(params.vertex_expanding_rule))};
  }
  return ExpandBoostedResult{.expanded_boosted = std::move(expanded_boosted)};
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
  return mongoose_match_s(graph, {}, params);
}

auto mongoose_match_s(const AdjacencyList<edge_probability_t>& graph, std::span<const vertex_id_t> seeds,
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

auto coarsen_wim_graph_by_match(const AdjacencyList<WIMEdge>& graph, const InvAdjacencyList<WIMEdge>& inv_graph,
                                std::span<const vertex_weight_t> vertex_weights, vertex_id_t n_groups,
                                std::span<const vertex_id_t> group_id, const CoarseningParams& params) noexcept
    -> rfl::Result<WIMCoarsenGraphBriefResult> {
  return coarsen_wim_graph_by_match_generic<WIMCoarsenGraphBriefResult>( //
      graph, inv_graph, vertex_weights, n_groups, group_id, params);
}

auto coarsen_wim_graph_by_match_d(const AdjacencyList<WIMEdge>& graph, const InvAdjacencyList<WIMEdge>& inv_graph,
                                  std::span<const vertex_weight_t> vertex_weights, vertex_id_t n_groups,
                                  std::span<const vertex_id_t> group_id, const CoarseningParams& params) noexcept
    -> rfl::Result<WIMCoarsenGraphDetailedResult> {
  return coarsen_wim_graph_by_match_generic<WIMCoarsenGraphDetailedResult>( //
      graph, inv_graph, vertex_weights, n_groups, group_id, params);
}

auto coarsen_wbim_graph_by_match(const AdjacencyList<WBIMEdge>& graph, const InvAdjacencyList<WBIMEdge>& inv_graph,
                                 std::span<const vertex_weight_t> vertex_weights, const VertexSet& seeds,
                                 vertex_id_t n_groups, std::span<const vertex_id_t> group_id,
                                 const CoarseningParams& params) noexcept -> rfl::Result<WBIMCoarsenGraphBriefResult> {
  return coarsen_wbim_graph_by_match_generic<WBIMCoarsenGraphBriefResult>( //
      graph, inv_graph, vertex_weights, seeds, n_groups, group_id, params);
}

auto coarsen_wbim_graph_by_match_d(const AdjacencyList<WBIMEdge>& graph, const InvAdjacencyList<WBIMEdge>& inv_graph,
                                   std::span<const vertex_weight_t> vertex_weights, const VertexSet& seeds,
                                   vertex_id_t n_groups, std::span<const vertex_id_t> group_id,
                                   const CoarseningParams& params) noexcept
    -> rfl::Result<WBIMCoarsenGraphDetailedResult> {
  return coarsen_wbim_graph_by_match_generic<WBIMCoarsenGraphDetailedResult>( //
      graph, inv_graph, vertex_weights, seeds, n_groups, group_id, params);
}

auto coarsen_wim_graph(const AdjacencyList<WIMEdge>& graph, const InvAdjacencyList<WIMEdge>& inv_graph,
                       std::span<const vertex_weight_t> vertex_weights, const CoarseningParams& params) noexcept
    -> rfl::Result<WIMCoarsenGraphBriefResult> {
  return coarsen_wim_graph_framework<WIMCoarsenGraphBriefResult>(graph, inv_graph, vertex_weights, params);
}

auto coarsen_wim_graph_d(const AdjacencyList<WIMEdge>& graph, const InvAdjacencyList<WIMEdge>& inv_graph,
                         std::span<const vertex_weight_t> vertex_weights, const CoarseningParams& params) noexcept
    -> rfl::Result<WIMCoarsenGraphDetailedResult> {
  return coarsen_wim_graph_framework<WIMCoarsenGraphDetailedResult>(graph, inv_graph, vertex_weights, params);
}

auto coarsen_wbim_graph(const AdjacencyList<WBIMEdge>& graph, const InvAdjacencyList<WBIMEdge>& inv_graph,
                        std::span<const vertex_weight_t> vertex_weights, const VertexSet& seeds,
                        const CoarseningParams& params) noexcept -> rfl::Result<WBIMCoarsenGraphBriefResult> {
  return coarsen_wbim_graph_framework<WBIMCoarsenGraphBriefResult>( //
      graph, inv_graph, vertex_weights, seeds, params);
}

auto coarsen_wbim_graph_d(const AdjacencyList<WBIMEdge>& graph, const InvAdjacencyList<WBIMEdge>& inv_graph,
                          std::span<const vertex_weight_t> vertex_weights, const VertexSet& seeds,
                          const CoarseningParams& params) noexcept -> rfl::Result<WBIMCoarsenGraphDetailedResult> {
  return coarsen_wbim_graph_framework<WBIMCoarsenGraphDetailedResult>( //
      graph, inv_graph, vertex_weights, seeds, params);
}

// ---- Step 4 ----

auto select_best_seed_in_group(const WIMCoarsenedVertexDetails& v) noexcept -> SelectBestSeedResult {
  auto m = v.n_members();
  if (m == 1) {
    return {.index_in_group = 0};
  }
  BOOST_ASSERT(m <= 3);
  auto& p_out = v.heuristics_out;
  auto& p_out_seed = v.heuristics_out_seed;
  auto to_estimated_gain = [&](size_t j0) {
    if (m == 2) {
      auto j1 = get_next_j<2>(j0);
      return p_out_seed[j0] + v.p_seed_internal[j0][j1] * p_out[j1];
    }
    BOOST_ASSERT(m == 3);
    auto [j1, j2] = get_next_j<3>(j0);
    auto p1 = at_least_1_probability(                      //
        v.p_seed_internal[j0][j1],                         // j0 -> j1
        v.p_seed_internal[j0][j2] * v.p_internal[j2][j1]); // j0 -> j2 -> j1
    auto p2 = at_least_1_probability(                      //
        v.p_seed_internal[j0][j2],                         // j0 -> j2
        v.p_seed_internal[j0][j1] * v.p_internal[j1][j2]); // j0 -> j1 -> j2
    return p_out_seed[j0] + p1 * p_out[j1] + p2 * p_out[j2];
  };
  auto estimated_gain = [&]() -> StaticVector<edge_probability_t, MAX_N_MEMBERS> {
    auto view = range(m) | views::transform(to_estimated_gain);
    return {view.begin(), view.end()};
  }();
  MYLOG_FMT_TRACE("Best seed selection: estimated gain for group {} = {::.4f}", v.members, estimated_gain);
  return {.index_in_group = static_cast<size_t>(ranges::max_element(estimated_gain) - estimated_gain.begin())};
}

namespace {
auto select_best_boosted_in_group_as_source(const WBIMCoarsenedVertexDetails& v) noexcept -> SelectBestBoostedResult {
  auto m = v.n_members();
  if (m == 1) {
    return {.index_in_group = 0};
  }
  BOOST_ASSERT(m <= 3);
  auto to_estimated_gain = [&](size_t j0) -> edge_probability_t {
    auto delta_h = v.heuristics_in_boost[j0] - v.heuristics_in[j0];
    if (m == 2) {
      auto j1 = get_next_j<2>(j0);
      return delta_h * (v.heuristics_out[j0] + v.p_internal[j0][j1] * v.heuristics_out[j1]);
    }
    BOOST_ASSERT(m == 3);
    auto [j1, j2] = get_next_j<3>(j0);
    auto p1 = at_least_1_probability(v.p_internal[j0][j1], v.p_internal[j0][j2] * v.p_internal[j2][j1]);
    auto p2 = at_least_1_probability(v.p_internal[j0][j2], v.p_internal[j0][j1] * v.p_internal[j1][j2]);
    return delta_h * (v.heuristics_out[j0] + p1 * v.heuristics_out[j1] + p2 * v.heuristics_out[j2]);
  };
  auto estimated_gain = [&]() -> StaticVector<edge_probability_t, MAX_N_MEMBERS> {
    auto view = range(m) | views::transform(to_estimated_gain);
    return {view.begin(), view.end()};
  }();
  MYLOG_FMT_TRACE("Best boosted vertex selection: estimated gain for group {} = {::.4f}", v.members, estimated_gain);
  return {.index_in_group = static_cast<size_t>(ranges::max_element(estimated_gain) - estimated_gain.begin())};
}

auto select_best_boosted_in_group_as_target(const WBIMCoarsenedVertexDetails& v) noexcept -> SelectBestBoostedResult {
  auto m = v.n_members();
  if (m == 1) {
    return {.index_in_group = 0};
  }
  BOOST_ASSERT(m <= 3);
  auto to_estimated_gain = [&](size_t j0) -> edge_probability_t {
    auto delta_h = v.heuristics_in_boost[j0] - v.heuristics_in[j0];
    if (m == 2) {
      auto j1 = get_next_j<2>(j0);
      auto delta_p_1 = v.p_boost_internal[j1][j0] - v.p_internal[j1][j0];
      // (1) InNeighbor(j0) -> [j0], (2) InNeighbor(j1) -> ... -> [j0]
      auto delta_in = delta_h + v.heuristics_in[j1] * delta_p_1;
      return delta_in * v.heuristics_out[j0];
    }
    BOOST_ASSERT(m == 3);
    auto [j1, j2] = get_next_j<3>(j0);
    // After - Before j0 is boosted: (1) j1 -> [j0]; (2) j1 -> j2 -> [j0]
    auto delta_p_1 =
        at_least_1_probability(v.p_boost_internal[j1][j0], v.p_internal[j1][j2] * v.p_boost_internal[j2][j0]) -
        at_least_1_probability(v.p_internal[j1][j0], v.p_internal[j1][j2] * v.p_internal[j2][j0]);
    // After - Before j0 is boosted: (1) j3 -> [j0]; (2) j2 -> j1 -> [j0]
    auto delta_p_2 =
        at_least_1_probability(v.p_boost_internal[j2][j0], v.p_internal[j2][j1] * v.p_boost_internal[j1][j0]) -
        at_least_1_probability(v.p_internal[j2][j0], v.p_internal[j2][j1] * v.p_internal[j1][j0]);
    // (1) InNeighbor(j0) -> [j0], (2) InNeighbor(j1) -> ... -> [j0], (3) InNeighbor(j2) -> ... -> [j0]
    auto delta_in = delta_h + v.heuristics_in[j1] * delta_p_1 + v.heuristics_in[j2] * delta_p_2;
    return delta_in * v.heuristics_out[j0];
  };
  auto estimated_gain = [&]() -> StaticVector<edge_probability_t, MAX_N_MEMBERS> {
    auto view = range(m) | views::transform(to_estimated_gain);
    return {view.begin(), view.end()};
  }();
  MYLOG_FMT_TRACE("Best boosted vertex selection: estimated gain for group {} = {::.4f}", v.members, estimated_gain);
  return {.index_in_group = static_cast<size_t>(ranges::max_element(estimated_gain) - estimated_gain.begin())};
}
} // namespace

auto select_best_boosted_in_group(const WBIMCoarsenedVertexDetails& v, BoostedSelectionRule rule) noexcept
    -> SelectBestBoostedResult {
  if (rule == BoostedSelectionRule::AS_SOURCE) {
    return select_best_boosted_in_group_as_source(v);
  } else {
    BOOST_ASSERT(rule == BoostedSelectionRule::AS_TARGET);
    return select_best_boosted_in_group_as_target(v);
  }
}

auto expand_wim_seed_vertices(const AdjacencyList<WIMEdge>& graph, std::span<const vertex_weight_t> vertex_weights,
                              const WIMCoarsenGraphBriefResult& coarsening_result,
                              std::span<const vertex_id_t> coarsened_seeds, const ExpandingParams& params) noexcept
    -> rfl::Result<ExpandSeedResult> {
  return expand_wim_seed_vertices_impl(graph, vertex_weights, coarsening_result, coarsened_seeds, params);
}

auto expand_wim_seed_vertices_d(const AdjacencyList<WIMEdge>& graph, std::span<const vertex_weight_t> vertex_weights,
                                const WIMCoarsenGraphDetailedResult& coarsening_result,
                                std::span<const vertex_id_t> coarsened_seeds, const ExpandingParams& params) noexcept
    -> rfl::Result<ExpandSeedResult> {
  return expand_wim_seed_vertices_impl(graph, vertex_weights, coarsening_result, coarsened_seeds, params);
}

auto expand_wbim_boosted_vertices(const AdjacencyList<WBIMEdge>& graph, std::span<const vertex_weight_t> vertex_weights,
                                  const WBIMCoarsenGraphBriefResult& coarsening_result,
                                  std::span<const vertex_id_t> coarsened_boosted,
                                  const ExpandingParams& params) noexcept -> rfl::Result<ExpandBoostedResult> {
  // Note: graph and vertex_weights are unused yet, reserved for further implementation.
  return expand_wbim_boosted_vertices_impl(coarsening_result, coarsened_boosted, params);
}

auto expand_wbim_boosted_vertices_d(const AdjacencyList<WBIMEdge>& graph,
                                    std::span<const vertex_weight_t> vertex_weights,
                                    const WBIMCoarsenGraphDetailedResult& coarsening_result,
                                    std::span<const vertex_id_t> coarsened_boosted,
                                    const ExpandingParams& params) noexcept -> rfl::Result<ExpandBoostedResult> {
  // Note: graph and vertex_weights are unused yet, reserved for further implementation.
  return expand_wbim_boosted_vertices_impl(coarsening_result, coarsened_boosted, params);
}
