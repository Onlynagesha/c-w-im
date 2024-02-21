#pragma once

#include "graph_types.h"

inline auto make_sample_wim_graph() {
  auto edge_list = DirectedEdgeList<WIMEdge>{};
  edge_list.open_for_push_back();

  // (source, target, p)
  using EdgeTuple = std::tuple<vertex_id_t, vertex_id_t, edge_probability_t>;
  auto edges = std::vector<EdgeTuple>{
      {0, 1, 0.70}, {1, 2, 0.71}, {3, 4, 0.72}, {4, 5, 0.73}, {6, 7, 0.74}, {7, 8, 0.75}, {0, 3, 0.76}, {3, 6, 0.77},
      {1, 4, 0.78}, {4, 7, 0.79}, {2, 5, 0.80}, {5, 8, 0.81}, {4, 1, 0.30}, {7, 4, 0.31}, {4, 3, 0.32}, {5, 4, 0.33}};
  for (auto [u, v, p] : edges) {
    edge_list.push_back(u, v, {.p = p, .p_seed = p + 0.1_ep});
  }
  edge_list.close_for_push_back();
  return std::tuple{AdjacencyList<WIMEdge>{edge_list}, InvAdjacencyList<WIMEdge>{edge_list}};
}
