#pragma once

#include "graph_types.h"
#include "utils/boost_assert.h"
#include "utils/dynamic_bitset.h"
#include "utils/graph.h"
#include <rfl/Size.hpp>
#include <rfl/Validator.hpp>
#include <rfl/comparisons.hpp>
#include <ylt/easylog/record.hpp>

template <class E>
struct AdjacencyListPair {
  AdjacencyList<E> adj_list;
  InvAdjacencyList<E> inv_adj_list;
  std::vector<vertex_weight_t> vertex_weights;

  auto in_degree(vertex_id_t v) const {
    BOOST_ASSERT_MSG(v >= 0 && v < graph::num_vertices(inv_adj_list), "v is out of range [0, n).");
    return graph::degree(inv_adj_list, v);
  }

  auto out_degree(vertex_id_t v) const {
    BOOST_ASSERT_MSG(v >= 0 && v < graph::num_vertices(inv_adj_list), "v is out of range [0, n).");
    return graph::degree(adj_list, v);
  }

  auto degree(vertex_id_t v) const {
    return in_degree(v) + out_degree(v);
  }
};

struct VertexSet {
  std::vector<vertex_id_t> vertex_list;
  DynamicBitset mask;

  explicit VertexSet(vertex_id_t n) : vertex_list(), mask(n) {}

  VertexSet(vertex_id_t n, std::span<const vertex_id_t> vertices)
      : vertex_list(vertices.begin(), vertices.end()), mask(n) {
    for (auto v : vertices) {
      BOOST_ASSERT_MSG(v >= 0 && v < n, "Vertex index out of range [0, n)");
      mask.set(v);
    }
  }

  VertexSet(vertex_id_t n, std::initializer_list<vertex_id_t> vertices)
      : VertexSet(n, {vertices.begin(), vertices.end()}) {}

  auto num_vertices_in_whole_graph() const -> vertex_id_t {
    return static_cast<vertex_id_t>(mask.size());
  }

  auto size() const -> vertex_id_t {
    return static_cast<vertex_id_t>(vertex_list.size());
  }

  auto contains(vertex_id_t v) const -> bool {
    BOOST_ASSERT_MSG(v >= 0 && v < num_vertices_in_whole_graph(), "Vertex index out of range [0, n)");
    return mask.test(v);
  }
};

struct RRSketchSet {
  const InvAdjacencyList<WIMEdge>* inv_graph;
  std::discrete_distribution<vertex_id_t> center_distribution;
  // sketches[i] = { v1, v2... }, all the vertices in the i-th RR-sketch
  std::vector<std::vector<vertex_id_t>> sketches;
  // inv_sketches[v] = { i1, i2 ... }, all the RR-sketches that contain v
  std::vector<std::vector<size_t>> inv_sketches;

  explicit RRSketchSet(const InvAdjacencyList<WIMEdge>* inv_graph, std::span<const vertex_weight_t> vertex_weights)
      : inv_graph(inv_graph), center_distribution(vertex_weights.begin(), vertex_weights.end()), sketches(),
        inv_sketches() {
    BOOST_ASSERT_MSG(vertex_weights.size() == graph::num_vertices(*inv_graph),
                     "Mismatch between # of vertices in the graph and # of vertices in the weight list.");
    auto n = graph::num_vertices(*inv_graph);
    inv_sketches.assign(n, std::vector<size_t>{});
  }

  auto num_vertices() const -> vertex_id_t {
    return static_cast<vertex_id_t>(inv_sketches.size());
  }

  auto num_sketches() const -> size_t {
    return sketches.size();
  }

  auto sketch_sizes() const {
    return sketches | views::transform(ranges::size);
  }

  auto average_sketch_size() const -> double {
    BOOST_ASSERT_MSG(!sketches.empty(), "Requires at least 1 RR-sketch to exist on calculating average size.");
    return 1.0 * accumulate_sum(sketch_sizes()) / sketches.size();
  }

  auto ratio_of_single_vertex_sketch() const -> double {
    BOOST_ASSERT_MSG(!sketches.empty(), "Requires at least 1 RR-sketch to exist on calculating average size.");
    return 1.0 * ranges::count(sketch_sizes(), 1) / sketches.size();
  }

  auto percentage_of_single_vertex_sketch() const -> double {
    return 100.0 * ratio_of_single_vertex_sketch();
  }

  auto rr_sketch_total_size_bytes() const -> size_t {
    auto res = sketches.capacity() * sizeof(decltype(sketches)::value_type) +
               inv_sketches.capacity() * sizeof(decltype(inv_sketches)::value_type);
    for (const auto& s : sketches) {
      res += s.capacity() * sizeof(vertex_id_t);
    }
    for (const auto& is : inv_sketches) {
      res += is.capacity() * sizeof(size_t);
    }
    return res;
  }

  auto rr_sketch_total_size_str() const -> std::string {
    constexpr auto UNITS = std::array{"Bytes", "KibiBytes", "Mebibytes", "Gibibytes"};
    auto value = static_cast<double>(rr_sketch_total_size_bytes());
    auto unit_index = 0;
    while (value >= 1024.0 && unit_index + 1 < UNITS.size()) {
      value /= 1024.0;
      unit_index += 1;
    }
    return fmt::format("{:.3f} {}", value, UNITS[unit_index]);
  }

  // Appends r new RR-sketches
  auto append(size_t n_sketches) noexcept -> void;

  // Selects k seed vertices
  auto select(vertex_id_t k) noexcept -> std::vector<vertex_id_t>;

  auto dump() noexcept -> std::string;

private:
  // Appends a new RR-sketch
  auto append_single(std::span<vertex_id_t> vertices) noexcept -> void;
};

auto read_wim_graph_data(const std::string& input_file) noexcept -> rfl::Result<AdjacencyListPair<WIMEdge>>;

auto read_wbim_graph_data(const std::string& input_file) noexcept -> rfl::Result<AdjacencyListPair<WBIMEdge>>;

auto wim_simulate(const AdjacencyList<WIMEdge>& graph, const VertexSet& seeds, uint64_t try_count) noexcept
    -> rfl::Result<double>;

auto wim_simulate_w(const AdjacencyList<WIMEdge>& graph, std::span<const vertex_weight_t> vertex_weights,
                    const VertexSet& seeds, uint64_t try_count) noexcept -> rfl::Result<double>;

auto wbim_simulate(const AdjacencyList<WBIMEdge>& graph, const VertexSet& seeds, const VertexSet& boosted_vertices,
                   uint64_t try_count) noexcept -> rfl::Result<double>;

auto wbim_simulate_w(const AdjacencyList<WBIMEdge>& graph, std::span<const vertex_weight_t> vertex_weights,
                     const VertexSet& seeds, const VertexSet& boosted_vertices, uint64_t try_count) noexcept
    -> rfl::Result<double>;

inline auto wim_simulate_s(const AdjacencyList<WIMEdge>& graph, std::span<const vertex_id_t> seeds, uint64_t try_count)
    -> rfl::Result<double> {
  auto seed_set = VertexSet{graph::num_vertices(graph), seeds};
  return wim_simulate(graph, seed_set, try_count);
}

inline auto wim_simulate_w_s(const AdjacencyList<WIMEdge>& graph, std::span<const vertex_weight_t> vertex_weights,
                             std::span<const vertex_id_t> seeds, uint64_t try_count) noexcept -> rfl::Result<double> {
  auto seed_set = VertexSet{graph::num_vertices(graph), seeds};
  return wim_simulate_w(graph, vertex_weights, seed_set, try_count);
}

inline auto wbim_simulate_s(const AdjacencyList<WBIMEdge>& graph, std::span<const vertex_id_t> seeds,
                            std::span<const vertex_id_t> boosted_vertices, uint64_t try_count) -> rfl::Result<double> {
  auto n = graph::num_vertices(graph);
  auto seed_set = VertexSet{n, seeds};
  auto boosted_set = VertexSet{n, boosted_vertices};
  return wbim_simulate(graph, seed_set, boosted_set, try_count);
}

inline auto wbim_simulate_w_s(const AdjacencyList<WBIMEdge>& graph, std::span<const vertex_weight_t> vertex_weights,
                              std::span<const vertex_id_t> seeds, std::span<const vertex_id_t> boosted_vertices,
                              uint64_t try_count) -> rfl::Result<double> {
  auto n = graph::num_vertices(graph);
  auto seed_set = VertexSet{n, seeds};
  auto boosted_set = VertexSet{n, boosted_vertices};
  return wbim_simulate_w(graph, vertex_weights, seed_set, boosted_set, try_count);
}
