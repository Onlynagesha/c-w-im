#pragma once

#include "graph_types.h"

namespace sample_graph_details {
using EdgeTuple = std::tuple<vertex_id_t, vertex_id_t, edge_probability_t>;

// (source, target, p)
inline const auto SAMPLE_EDGES_1 = std::vector<EdgeTuple>{
    {0, 1, 0.30}, {1, 2, 0.31}, {3, 4, 0.32}, {4, 5, 0.33}, {6, 7, 0.34}, {7, 8, 0.35}, {0, 3, 0.41},
    {3, 6, 0.44}, {1, 4, 0.42}, {4, 7, 0.45}, {2, 5, 0.43}, {5, 8, 0.46}, {4, 1, 0.20}, {7, 4, 0.21},
    {4, 3, 0.22}, {5, 4, 0.23}, {9, 0, 0.4},  {9, 3, 0.3},  {9, 6, 0.2},  {2, 0, 0.10}, {8, 6, 0.11}};

// (source, target, p)
inline const auto SAMPLE_EDGES_1_LEFT =
    std::vector<EdgeTuple>{{0, 1, 0.30}, {2, 3, 0.32},       {4, 5, 0.34},       {0, 2, 0.41},
                           {2, 4, 0.44}, {6, 0, 0.45045045}, {6, 2, 0.27027027}, {6, 4, 0.27927927}};

// (source, target, p)
inline const auto SAMPLE_EDGES_1_RIGHT =
    std::vector<EdgeTuple>{{0, 1, 0.31}, {2, 3, 0.33}, {4, 5, 0.35},   {0, 2, 0.42},     {2, 4, 0.45},
                           {2, 0, 0.20}, {4, 2, 0.21}, {6, 0, 0.3125}, {6, 2, 0.333333}, {6, 4, 0.354167}};

inline auto SAMPLE_EDGES_2 =
    std::vector<EdgeTuple>{{1, 2, 0.3}, {1, 3, 0.4}, {3, 2, 0.5}, {2, 4, 0.6}, {3, 5, 0.7}, {4, 5, 0.8}};

inline auto SAMPLE_EDGES_3 = std::vector<EdgeTuple>{
    {1, 0, 0.20},   {2, 0, 0.21},   {3, 0, 0.22},   {4, 2, 0.25},   {5, 2, 0.26},  {6, 2, 0.27},
    {7, 3, 0.28},   {9, 3, 0.29},   {8, 5, 0.30},   {10, 9, 0.31},  {13, 9, 0.32}, {14, 13, 0.35},
    {11, 13, 0.36}, {12, 11, 0.37}, {15, 13, 0.40}, {16, 15, 0.41}, {17, 1, 0.90}};

inline auto make_sample_wim_graph(std::span<const EdgeTuple> edges) {
  auto edge_list = DirectedEdgeList<WIMEdge>{};
  edge_list.open_for_push_back();
  for (auto [u, v, p] : edges) {
    edge_list.push_back(u, v, {.p = p, .p_seed = std::min(2.0_ep * p, 1.0_ep)});
  }
  edge_list.close_for_push_back();
  return std::tuple{AdjacencyList<WIMEdge>{edge_list}, InvAdjacencyList<WIMEdge>{edge_list}};
}

inline auto make_sample_wbim_graph(std::span<const EdgeTuple> edges) {
  auto edge_list = DirectedEdgeList<WBIMEdge>{};
  edge_list.open_for_push_back();
  for (auto [u, v, p] : edges) {
    edge_list.push_back(u, v, {.p = p, .p_boost = std::min(2.0_ep * p, 1.0_ep)});
  }
  edge_list.close_for_push_back();
  return std::tuple{AdjacencyList<WBIMEdge>{edge_list}, InvAdjacencyList<WBIMEdge>{edge_list}};
}
} // namespace sample_graph_details

inline auto make_sample_wim_graph_1() {
  return sample_graph_details::make_sample_wim_graph(sample_graph_details::SAMPLE_EDGES_1);
}

inline auto make_sample_wbim_graph_1() {
  return sample_graph_details::make_sample_wbim_graph(sample_graph_details::SAMPLE_EDGES_1);
}

inline auto make_sample_wim_graph_1_left() {
  return sample_graph_details::make_sample_wim_graph(sample_graph_details::SAMPLE_EDGES_1_LEFT);
}

inline auto make_sample_wim_graph_1_right() {
  return sample_graph_details::make_sample_wim_graph(sample_graph_details::SAMPLE_EDGES_1_RIGHT);
}

inline auto make_sample_wim_graph_2() {
  return sample_graph_details::make_sample_wim_graph(sample_graph_details::SAMPLE_EDGES_2);
}

inline auto make_sample_wbim_graph_3() {
  return sample_graph_details::make_sample_wbim_graph(sample_graph_details::SAMPLE_EDGES_3);
}
