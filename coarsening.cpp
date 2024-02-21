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
    return std::numeric_limits<edge_probability_t>::lowest();
  case NeighborMatchRule::LEM_P_MAX:
  case NeighborMatchRule::LEM_P_PRODUCT:
    return std::numeric_limits<edge_probability_t>::max();
  default:
    BOOST_ASSERT(!"Unreachable branch.");
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

template <is_edge_property E>
auto merge_edge_to_undirected_generic(AdjacencyList<E>& graph, const CoarseningParams& params)
    -> AdjacencyList<edge_probability_t> {
  // Ensures |V| is unchanged
  auto edge_list = UndirectedEdgeList<edge_probability_t>{graph::num_vertices(graph)};
  edge_list.open_for_push_back();

  using VertexPair = std::pair<vertex_id_t, vertex_id_t>;
  auto edges_map = std::map<VertexPair, edge_probability_t>{};
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

// The probability of spread success from gu.members[j0] to any member of gv
template <size_t N, bool Seed>
auto get_p_across_groups(size_t j0, const CoarsenedEdgeDetails& edge, const CoarsenedVertexDetails& gu,
                         const CoarsenedVertexDetails& gv, const CoarseningParams& params) -> edge_probability_t {
  auto [p_cross_first, p_internal_first] =
      Seed ? std::tuple{&edge.p_seed_cross, &gu.p_seed_internal} : std::tuple{&edge.p_cross, &gu.p_internal};
  if constexpr (N == 2) {
    auto paths = StaticVector<edge_probability_t, 6>{};
    auto j1 = (j0 + 1) % 2;
    for (auto i : range(gv.members.size())) {
      paths.push_back((*p_cross_first)[j0][i]);                           // (1) j0 -> i
      paths.push_back((*p_internal_first)[j0][j1] * edge.p_cross[j1][i]); // (2) j0 -> j1 -> i
    }
    return at_least_1_probability_of_range(paths);
  } else if constexpr (N == 3) {
    auto paths = StaticVector<edge_probability_t, 15>{};
    auto j1 = (j0 + 1) % 3;
    auto j2 = (j0 + 2) % 3;
    for (auto i : range(gv.members.size())) {
      paths.push_back((*p_cross_first)[j0][i]);                           // (1) j0 -> i
      paths.push_back((*p_internal_first)[j0][j1] * edge.p_cross[j1][i]); // (2) j0 -> j1 -> i
      paths.push_back((*p_internal_first)[j0][j2] * edge.p_cross[j2][i]); // (3) j0 -> j2 -> i
      // (4,5) j0 -> j1|j2 -> j2|j1 -> i
      paths.push_back((*p_internal_first)[j0][j1] * gu.p_internal[j1][j2] * edge.p_cross[j2][i]);
      paths.push_back((*p_internal_first)[j0][j2] * gu.p_internal[j2][j1] * edge.p_cross[j1][i]);
    }
    return at_least_1_probability_of_range(paths);
  } else {
    static_assert(rfl::always_false_v<std::make_index_sequence<N>>, "Invalid group size.");
  }
}

template <size_t N>
  requires(N == 2 || N == 3)
auto merge_coarsened_wim_graph_edge_common(CoarsenedEdgeDetails& dest, const CoarsenedVertexDetails& gu,
                                           const CoarsenedVertexDetails& gv, const CoarseningParams& params) -> void {
  BOOST_ASSERT(gu.members.size() == N);
  // Part (1): merged.p
  auto heuristics_sum = accumulate_sum(gu.heuristics_in);
  if (heuristics_sum <= 0.0) {
    // Cannot serve as non-seed coarsened vertex since the group is isolated
    // (with no in-edge from other coarsened vertices)
    dest.merged.p = 0.0;
  } else {
    auto p_paths = std::array<edge_probability_t, N>{};
    for (auto j0 : range(N)) {
      p_paths[j0] = get_p_across_groups<N, false>(j0, dest, gu, gv, params);
    }
    ELOGFMT(TRACE, "p_paths from group {} to group {} = {::.4f}", gu.members, gv.members, p_paths);
    ELOGFMT(TRACE, "normalized in-heuristic of group {} = {::.4f}", gu.members,
            gu.heuristics_in | views::transform(LAMBDA_1(_1 / heuristics_sum)));

    auto to_combined = views::transform(LAMBDA_1(gu.heuristics_in[_1] * p_paths[_1] / heuristics_sum));
    if (params.group_path_probability_rule == GroupPathProbabilityRule::P_MERGED) {
      auto to_probability_in = views::transform(LAMBDA_1(gu.heuristics_in[_1] / heuristics_sum));
      dest.merged.p = at_least_1_probability_of_range(range(N) | to_combined) /
                      at_least_1_probability_of_range(range(N) | to_probability_in);
    } else {
      dest.merged.p = accumulate_sum(range(N) | to_combined);
    }
  }
  // Part (2): merged.p_boost, which only considers the best candidate inside the group gu
  auto estimated_gain_as_seed = std::array<edge_probability_t, N>{};
  for (auto j0 : range(N)) {
    estimated_gain_as_seed[j0] = get_p_across_groups<N, true>(j0, dest, gu, gv, params);
  }
  auto best_at = ranges::max_element(estimated_gain_as_seed);
  ELOGFMT(TRACE, "Estimated gain as seed of group {}: {} (best vertex as seed: {})", gu.members, estimated_gain_as_seed,
          gu.members[best_at - estimated_gain_as_seed.begin()]);
  dest.merged.p_seed = *best_at;
}

auto merge_coarsened_wim_graph_edge(CoarsenedEdgeDetails& dest, const CoarsenedVertexDetails& gu,
                                    const CoarsenedVertexDetails& gv, const CoarseningParams& params) -> void {
  switch (gu.members.size()) {
  case 1:
    // At least one of u -> v0, u -> v1, u -> v2
    dest.merged.p = at_least_1_probability_of_range(dest.p_cross[0] | views::take(gv.members.size()));
    dest.merged.p_seed = at_least_1_probability_of_range(dest.p_seed_cross[0] | views::take(gv.members.size()));
    break;
  case 2:
    merge_coarsened_wim_graph_edge_common<2>(dest, gu, gv, params);
    break;
  case 3:
    merge_coarsened_wim_graph_edge_common<3>(dest, gu, gv, params);
    break;
  default:
    BOOST_ASSERT_MSG(false, "Unreachable branch.");
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
  for (auto [u, v] : graph::make_edge_range<>(graph)) {
    auto [gu, gv] = details.to_group_id(u, v);
    if (gu == gv) {
      continue; // Skips the edges (u, v) whose two ends locate in the same group
    }
    auto it = res.find({gu, gv});
    BOOST_ASSERT(it != res.end());
    if (it->second.merged.p == -1.0_ep) {
      merge_coarsened_wim_graph_edge(it->second, details.groups[gu], details.groups[gv], params);
    }
  }
  return res;
}

auto get_coarsened_wim_graph_vertex_weights(const CoarseningDetails& details,
                                            std::span<const vertex_weight_t> vertex_weights,
                                            const CoarseningParams& params) -> std::vector<vertex_weight_t> {
  auto coarsened_vertex_weights = std::vector<vertex_weight_t>(details.n_groups, 0.0);
  // vc = Coarsened vertex index in range [0, Nc)
  for (auto [vc, gc] : details.groups | views::enumerate) {
    auto h_in_sum = accumulate_sum(gc.heuristics_in);
    auto& dest = coarsened_vertex_weights[vc];

    if (params.vertex_weight_rule == VertexWeightRule::SUM) {
      dest = accumulate_sum(gc.members | views::transform(LAMBDA_1(vertex_weights[_1])));
    } else if (h_in_sum <= 0.0) {
      // Special case for AVERAGE and AVERAGE_BY_PATHS: Takes arithmetic average if no in-heuristic
      dest = accumulate_sum(gc.members | views::transform(LAMBDA_1(vertex_weights[_1]))) / gc.members.size();
    } else if (params.vertex_weight_rule == VertexWeightRule::AVERAGE) {
      // AVERAGE: Weighted average by normalized in-heuristic
      auto to_vertex_weight = LAMBDA_1(vertex_weights[gc.members[_1]] * gc.heuristics_in[_1] / h_in_sum);
      dest = accumulate_sum(range(gc.members.size()) | views::transform(to_vertex_weight));
    } else if (gc.members.size() == 1) {
      // AVERAGE_BY_PATHS for group of size 1
      dest = vertex_weights[gc.members[0]];
    } else if (gc.members.size() == 2) {
      // AVERAGE_BY_PATHS for group of size 2
      for (auto [i1, v1] : gc.members | views::enumerate) {
        // (1) v1; (2) v0 -> v1
        auto i0 = (i1 + 1) % 2;
        auto paths = {gc.heuristics_in[i1] / h_in_sum, gc.heuristics_in[i0] * gc.p_internal[i0][i1] / h_in_sum};
        dest += vertex_weights[v1] * at_least_1_probability_of_range(paths);
      }
    } else {
      // AVERAGE_BY_PATHS for group of size 3
      BOOST_ASSERT(gc.members.size() == 3);
      for (auto [i2, v2] : gc.members | views::enumerate) {
        auto i0 = (i2 + 1) % 3;
        auto i1 = (i2 + 2) % 3;
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

auto coarsen_wim_graph_w_impl(AdjacencyList<WIMEdge>& graph, InvAdjacencyList<WIMEdge>& inv_graph,
                              std::span<const vertex_weight_t> vertex_weights, vertex_id_t n_groups,
                              std::span<const vertex_id_t> group_id, const CoarseningParams& params) noexcept
    -> CoarsenGraphResult {
  ELOGFMT(DEBUG, "Parameters: {:4}", params);
  auto n = graph::num_vertices(graph);
  // d_in: unit weight by default
  auto details = CoarseningDetails{.n = n,
                                   .n_groups = n_groups,
                                   .group_id = std::vector(group_id.begin(), group_id.end()),
                                   .index_in_group = std::vector<vertex_id_t>(n),
                                   .groups = std::vector<CoarsenedVertexDetails>(n_groups)};
  ELOGFMT(DEBUG, "n = {}, n_groups = {}", details.n, details.n_groups);
  // Step 1: groups[].members & groups[].vertex_weights & groups[].heuristics_in & index_in_group
  // Vertex v is placed in group g
  for (auto [v, g] : views::enumerate(group_id)) {
    auto& cur_group = details.groups[g];
    details.index_in_group[v] = static_cast<vertex_id_t>(cur_group.members.size());
    cur_group.members.push_back(v);
    cur_group.vertex_weights.push_back(vertex_weights[v]);
    // neighbor_edges: list of pairs (u, w)
    auto get_heuristics = [&](vertex_id_t v, auto&& neighbor_edges, auto which_p) -> edge_probability_t {
      auto res = 0.0_ep;
      for (auto [u, w] : neighbor_edges) {
        if (group_id[u] == group_id[v]) {
          continue;
        }
        if (params.group_in_out_rule == GroupInOutRule::P) {
          res += w.*which_p;
        } else if (params.group_in_out_rule == GroupInOutRule::W) {
          res += w.*which_p * vertex_weights[u];
        } else {
          res += 1.0_ep;
        }
      }
      return res;
    };
    auto heuristics_in = get_heuristics(v, inv_graph[v], &WIMEdge::p); // In-neighbors
    cur_group.heuristics_in.push_back(heuristics_in);
    auto heuristics_out = get_heuristics(v, graph[v], &WIMEdge::p); // Out-neighbors
    cur_group.heuristics_out.push_back(heuristics_out);
    auto heuristics_out_seed = get_heuristics(v, graph[v], &WIMEdge::p_seed);
    cur_group.heuristics_out_seed.push_back(heuristics_out_seed);
  }
  // Step 2: groups[].p_internal
  for (auto [u, v, w] : graph::make_edge_range<0>(graph)) {
    auto g = group_id[u];
    if (g == group_id[v]) {
      auto [ju, jv] = details.to_index_in_group(u, v);
      details.groups[g].p_internal[ju][jv] = w.p;
      details.groups[g].p_seed_internal[ju][jv] = w.p_seed;
    }
  }
  // Step 3: Generates coarsened edges
  auto edge_details = get_coarsened_wim_graph_edges(graph, details, params);
  // Step 4: Generates coarsened vertex weights
  auto coarsened_vertex_weights = get_coarsened_wim_graph_vertex_weights(details, vertex_weights, params);
  // Finally, builds the result from all the components
  auto coarsened_edge_list = DirectedEdgeList<WIMEdge>{};
  coarsened_edge_list.open_for_push_back();
  for (const auto& [p, w] : edge_details) {
    coarsened_edge_list.push_back(get<0>(p), get<1>(p), w.merged);
  }
  coarsened_edge_list.close_for_push_back();

  return CoarsenGraphResult{.coarsened_graph = AdjacencyList<WIMEdge>{coarsened_edge_list},
                            .coarsened_inv_graph = InvAdjacencyList<WIMEdge>{coarsened_edge_list},
                            .coarsened_vertex_weights = std::move(coarsened_vertex_weights),
                            .details = std::move(details)};
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
  if (params.seed_merge_rule == SeedMergeRule::S_SINGLE) {
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
  if (params.seed_merge_rule == SeedMergeRule::S_SINGLE) {
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

auto do_coarsen_wim_graph(const AdjacencyList<WIMEdge>& graph, const InvAdjacencyList<WIMEdge>& inv_graph,
                          vertex_id_t n_groups, std::span<const vertex_id_t> group_id,
                          const CoarseningParams& params) noexcept -> CoarsenGraphResult {
  auto vertex_weights = std::vector<vertex_weight_t>(graph::num_vertices(graph), 1.0);
  return do_coarsen_wim_graph_w(graph, inv_graph, vertex_weights, n_groups, group_id, params);
}

auto do_coarsen_wim_graph_w(const AdjacencyList<WIMEdge>& graph, const InvAdjacencyList<WIMEdge>& inv_graph,
                            std::span<const vertex_weight_t> vertex_weights, vertex_id_t n_groups,
                            std::span<const vertex_id_t> group_id, const CoarseningParams& params) noexcept
    -> CoarsenGraphResult {
  // Workaround: const reference triggers a bug in nwgraph. The graph is not changed actually.
  auto& graph_non_const = const_cast<AdjacencyList<WIMEdge>&>(graph);
  auto& inv_graph_non_const = const_cast<InvAdjacencyList<WIMEdge>&>(inv_graph);
  return coarsen_wim_graph_w_impl(graph_non_const, inv_graph_non_const, vertex_weights, n_groups, group_id, params);
}

auto coarsen_wim_graph(const AdjacencyList<WIMEdge>& graph, const InvAdjacencyList<WIMEdge>& inv_graph,
                       const CoarseningParams& params) noexcept -> CoarsenGraphResult {
  auto vertex_weights = std::vector<vertex_weight_t>(graph::num_vertices(graph), 1.0);
  return coarsen_wim_graph_w(graph, inv_graph, vertex_weights, params);
}

auto coarsen_wim_graph_w(const AdjacencyList<WIMEdge>& graph, const InvAdjacencyList<WIMEdge>& inv_graph,
                         std::span<const vertex_weight_t> vertex_weights, const CoarseningParams& params) noexcept
    -> CoarsenGraphResult {
  // Step 1:
  ELOGFMT(DEBUG, "|V|, |E| of graph = {}; |V|, |E| of inv_graph = {}; Size of vertex_weights = {}", //
          graph_n_m(graph), graph_n_m(inv_graph), vertex_weights.size());
  auto bidir_graph = merge_wim_edge_to_undirected(graph, params);
  ELOGFMT(DEBUG, "|V|, |E| of bidir_graph = {}", graph_n_m(bidir_graph));
  // Step 2
  auto [n_groups, group_id] = mongoose_match(bidir_graph, params);
  ELOGFMT(DEBUG, "n_groups = {}; Size of group_id = {}", n_groups, group_id.size());
  // Step 3
  return do_coarsen_wim_graph_w(graph, inv_graph, vertex_weights, n_groups, group_id, params);
}

auto further_coarsen_wim_graph(const CoarsenGraphResult& last_result, const CoarseningParams& params) noexcept
    -> CoarsenGraphResult {
  return coarsen_wim_graph_w(last_result.coarsened_graph, last_result.coarsened_inv_graph,
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
    if (v.members.size() == 2) {
      auto j1 = (j0 + 1) % 2;
      res = v.heuristics_out_seed[j0] + v.p_seed_internal[j0][j1] * v.heuristics_out[j1];
    } else {
      BOOST_ASSERT(v.members.size() == 3);
      auto j1 = (j0 + 1) % 3;
      auto j2 = (j0 + 2) % 3;
      auto p1 = at_least_1_probability(v.p_seed_internal[j0][j1], v.p_seed_internal[j0][j2] * v.p_internal[j2][j1]);
      auto p2 = at_least_1_probability(v.p_seed_internal[j0][j2], v.p_seed_internal[j0][j1] * v.p_internal[j1][j2]);
      res = v.heuristics_out_seed[j0] + p1 * v.heuristics_out[j1] + p2 * v.heuristics_out[j2];
    }
    ELOGFMT(TRACE, "Best seed selection: estimated gain for #{} (index = {}) = {:.4f}", v.members[j0], j0, res);
    return res;
  };
  auto estimated_gain = [&]() {
    auto view = range(m) | views::transform(to_estimated_gain);
    return std::vector(view.begin(), view.end());
  }();
  return {.index_in_group = static_cast<size_t>(ranges::max_element(estimated_gain) - estimated_gain.begin())};
}

auto expand_wim_seed_vertices(const AdjacencyList<WIMEdge>& graph, const CoarsenGraphResult& coarsening_result,
                              std::span<const vertex_id_t> coarsened_seeds, const ExpandingParams& params) noexcept
    -> rfl::Result<ExpandSeedResult> {
  auto vertex_weights = std::vector<vertex_weight_t>(graph::num_vertices(graph), 1.0);
  return expand_wim_seed_vertices_w(graph, vertex_weights, coarsening_result, coarsened_seeds, params);
}

auto expand_wim_seed_vertices_w(const AdjacencyList<WIMEdge>& graph, std::span<const vertex_weight_t> vertex_weights,
                                const CoarsenGraphResult& coarsening_result,
                                std::span<const vertex_id_t> coarsened_seeds, const ExpandingParams& params) noexcept
    -> rfl::Result<ExpandSeedResult> {
  auto n = coarsening_result.details.n;         // # of vertices BEFORE coarsening
  auto nc = coarsening_result.details.n_groups; // # of vertices AFTER coarsening
  BOOST_ASSERT_MSG(ranges::size(coarsened_seeds) <= nc, "Expects more than nc vertices as seeds");
  BOOST_ASSERT_MSG(ranges::max(coarsened_seeds) < nc, "Indices of seeds are out of range [0, nc).");

  auto expanded_seeds = make_reserved_vector<vertex_id_t>(coarsened_seeds.size());
  auto make_expanded_seeds_as_s_local = [&]() {
    for (auto sc : coarsened_seeds) {
      auto& group = coarsening_result.details.groups[sc];
      auto best_seed_index = select_best_seed_in_group(group).index_in_group;
      expanded_seeds.push_back(group.members[best_seed_index]);
    }
  };
  ELOGFMT(DEBUG, "coarsened_seeds = {}", coarsened_seeds);
  if (params.seed_expansion_rule == SeedExpansionRule::S_LOCAL) {
    // S_LOCAL: Only considers information inside current group
    make_expanded_seeds_as_s_local();
  } else if (params.seed_expansion_rule == SeedExpansionRule::S_SIMULATIVE) {
    // S_SIMULATIVE: Picks one-by-one separately by simulation result of single expanded seed vertex
    for (auto sc : coarsened_seeds) {
      auto& group = coarsening_result.details.groups[sc];
      auto estimated_gain = StaticVector<double, CoarsenedVertexDetails::MAX_NUM_MEMBERS>{};
      for (auto v : group.members) {
        auto sim_res = wim_simulate_w(graph, vertex_weights, {n, {v}}, **params.simulation_try_count);
        if (!sim_res) {
          return *sim_res.error();
        }
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
    BOOST_ASSERT(params.seed_expansion_rule == SeedExpansionRule::S_ITERATIVE);
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
