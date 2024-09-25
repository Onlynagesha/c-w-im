#include "contrast_algorithms.h"
#include "dump.h"
#include "dump_utils.h"
#include "graph_connectivity.h"
#include "radix_sort.h"
#include "wim.h"

namespace {
struct VertexInGroup {
  vertex_id_t v;
  vertex_id_t v_group_index;

  constexpr auto operator<=>(const VertexInGroup& rhs) const -> std::weak_ordering {
    return v_group_index <=> rhs.v_group_index;
  }
  constexpr auto operator==(const VertexInGroup& rhs) const {
    return v_group_index == rhs.v_group_index;
  }
};

struct SCCPair {
  vertex_id_t v;
  vertex_id_t before;
  vertex_id_t after;

  constexpr auto operator<=>(const SCCPair& rhs) const -> std::weak_ordering {
    if (before != rhs.before) {
      return before <=> rhs.before;
    }
    return after <=> rhs.after;
  }

  constexpr auto operator==(const SCCPair& rhs) const {
    return before == rhs.before && after == rhs.after;
  }
};
} // namespace

auto wim_r_robust_scc(const AdjacencyList<WIMEdge>& graph, std::span<const vertex_weight_t> vertex_weights,
                      const RRobustSCCParams& params) -> rfl::Result<std::vector<vertex_id_t>> {
  auto [n, m] = graph_n_m(graph);
  if (vertex_weights.size() != 0 && vertex_weights.size() != n) {
    constexpr auto msg_pattern = "Invalid size of vertex weight list: {} expected, {} actual.";
    return rfl::Error{fmt::format(msg_pattern, n, vertex_weights.size())};
  }

  auto group_index = std::vector<vertex_id_t>(n);
  auto scc_pairs = make_reserved_vector<SCCPair>(n);
  auto edge_filter_fn = [](vertex_id_t, vertex_id_t, const WIMEdge& e) { return rand_bool(e.p); };

  auto n_coarsened = 0_vid;
  // Step 1: Groups vertices by r-robust SCCs
  for (auto iteration : range(params.r)) {
    auto [n_sccs, component_index] = strongly_connected_components(graph, edge_filter_fn);
    scc_pairs.clear();
    for (auto v : range(n)) {
      scc_pairs.push_back({.v = v, .before = group_index[v], .after = component_index[v]});
    }
    radix_sort_struct<&SCCPair::before, &SCCPair::after>(scc_pairs, 0_vid, n - 1);

    n_coarsened = 0_vid;
    for (auto chunk : scc_pairs | CHUNK_BY_VIEW(_1 == _2)) {
      for (const auto& p : chunk) {
        group_index[p.v] = n_coarsened;
      }
      n_coarsened += 1;
    }
    MYLOG_FMT_DEBUG("After {} iterations: coarsened |V| = {}", iteration + 1, n_coarsened);
    MYLOG_TRACE(DUMP_INDEX_ARRAY(group_index));
  }
  // Step 2: Gets coarsened vertex weights
  auto coarsened_vertex_weights = std::vector<vertex_weight_t>(n_coarsened);
  if (vertex_weights.size() != 0) {
    for (auto [v, w] : vertex_weights | views::enumerate) {
      coarsened_vertex_weights[group_index[v]] += w;
    }
  } else {
    for (auto v : vertices(graph)) {
      coarsened_vertex_weights[group_index[v]] += 1.0_vw;
    }
  }
  // Step 3: Gets coarsened edge probabilities
  auto p_map = std::map<std::pair<vertex_id_t, vertex_id_t>, edge_probability_t>{};
  for (auto u : vertices(graph)) {
    auto gu = group_index[u];
    for (auto [v, w] : graph[u]) {
      auto gv = group_index[v];
      if (gu == gv) {
        continue;
      }
      auto [iter, is_new_pair] = p_map.emplace(std::pair{gu, gv}, 1.0_ep - w.p);
      if (!is_new_pair) {
        iter->second *= (1.0_ep - w.p);
      }
    }
  }
  MYLOG_FMT_DEBUG("coarsened |E| = {}", p_map.size());
  // Step 4: Builds the coarsened graph
  auto [coarsened_adj_list, coarsened_inv_adj_list] = [&]() {
    auto coarsened_edge_list = DirectedEdgeList<WIMEdge>{};
    coarsened_edge_list.open_for_push_back();
    for (auto [v_pair, one_minus_p] : p_map) {
      auto [u, v] = v_pair;
      auto p = 1.0_ep - one_minus_p;
      coarsened_edge_list.push_back(u, v, {.p = p, .p_seed = p}); // p_seed is unused
    }
    coarsened_edge_list.close_for_push_back();
    return std::pair{AdjacencyList<WIMEdge>{coarsened_edge_list}, InvAdjacencyList<WIMEdge>{coarsened_edge_list}};
  }();
  MYLOG_FMT_TRACE("Coarsened graph:\n{}", dump_utils::dump_graph(coarsened_adj_list, 4));

  // Step 5: Solves by reverse-reachable sketching
  auto sketch_set = RRSketchSet{&coarsened_inv_adj_list, coarsened_vertex_weights};
  sketch_set.append(params.n_sketches);
  auto coarsened_seeds = sketch_set.select(params.k);

  // Step 6: Expands coarsened seeds to the original graph
  auto vertices_in_group = make_reserved_vector<VertexInGroup>(n);
  for (auto v : vertices(graph)) {
    vertices_in_group.push_back({.v = v, .v_group_index = group_index[v]});
  }
  radix_sort_struct<&VertexInGroup::v_group_index>(vertices_in_group, 0_vid, n_coarsened - 1);

  auto expanded_seeds = make_reserved_vector<vertex_id_t>(params.k);
  for (auto coarsened_s : coarsened_seeds) {
    auto chunk = ranges::equal_range(vertices_in_group, VertexInGroup{.v_group_index = coarsened_s});
    expanded_seeds.push_back(rand_element(chunk).v); // Selects the expanded seed randomly
  }
  // Finally
  return expanded_seeds;
}
