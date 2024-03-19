#pragma once

#include "graph_types.h"
#include "utils/dynamic_bitset.h"

struct VertexSet {
  std::vector<vertex_id_t> vertex_list;
  DynamicBitset mask;

  explicit VertexSet(vertex_id_t n) : vertex_list(), mask(n) {}

  VertexSet(vertex_id_t n, std::vector<vertex_id_t> vertices) : vertex_list(std::move(vertices)), mask(n) {
    for (auto v : vertex_list) {
      BOOST_ASSERT_MSG(v >= 0 && v < n, "Vertex index out of range [0, n)");
      mask.set(v);
    }
  }

  VertexSet(vertex_id_t n, std::span<const vertex_id_t> vertices)
      : vertex_list(vertices.begin(), vertices.end()), mask(n) {
    for (auto v : vertices) {
      BOOST_ASSERT_MSG(v >= 0 && v < n, "Vertex index out of range [0, n)");
      mask.set(v);
    }
  }

  VertexSet(vertex_id_t n, std::initializer_list<vertex_id_t> vertices)
      : VertexSet(n, std::span{vertices.begin(), vertices.end()}) {}

  auto n_vertices_in_whole_graph() const -> vertex_id_t {
    return static_cast<vertex_id_t>(mask.size());
  }

  auto size() const -> vertex_id_t {
    return static_cast<vertex_id_t>(vertex_list.size());
  }

  auto contains(vertex_id_t v) const -> bool {
    BOOST_ASSERT_MSG(v >= 0 && v < n_vertices_in_whole_graph(), "Vertex index out of range [0, n)");
    return mask.test(v);
  }
};
