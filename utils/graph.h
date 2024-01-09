#pragma once

#include <nwgraph/adjacency.hpp>
#include <nwgraph/edge_list.hpp>
#include <nwgraph/graph_concepts.hpp>

namespace graph = nw::graph;

constexpr auto DIRECTED = graph::directedness::directed;
constexpr auto UNDIRECTED = graph::directedness::undirected;

using vertex_id_t = graph::default_vertex_id_type;
using edge_id_t = graph::default_index_t;

template <class... Attributes>
using DirectedEdgeList = graph::edge_list<DIRECTED, Attributes...>;

template <class... Attributes>
using UndirectedEdgeList = graph::edge_list<UNDIRECTED, Attributes...>;

template <class... Attributes>
using AdjacencyList = graph::adjacency<0, Attributes...>;

template <class... Attributes>
using InvAdjacencyList = graph::adjacency<1, Attributes...>;

inline constexpr auto operator""_vid(unsigned long long v) {
  return static_cast<vertex_id_t>(v);
}

inline constexpr auto operator""_eid(unsigned long long v) {
  return static_cast<edge_id_t>(v);
}

template <graph::graph Graph>
inline auto vertices(const Graph& g) {
  auto n = graph::num_vertices(g);
  return std::views::iota(static_cast<decltype(n)>(0), n);
}
