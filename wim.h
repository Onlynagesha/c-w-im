#pragma once

#include "graph_types.h"
#include "utils/boost_assert.h"
#include "utils/dynamic_bitset.h"
#include "utils/graph.h"
#include <memory_resource>
#include <rfl/Size.hpp>
#include <rfl/Validator.hpp>
#include <rfl/comparisons.hpp>
#include <ylt/easylog/record.hpp>

template <class E>
struct AdjacencyListPair {
  AdjacencyList<E> adj_list;
  InvAdjacencyList<E> inv_adj_list;
  std::vector<vertex_weight_t> vertex_weights;
};

struct WIMParams {
  std::string input_file;
  std::vector<size_t> num_sketches;
  std::vector<vertex_id_t> num_seeds;
  rfl::Validator<uint64_t, rfl::Minimum<1>> simulation_try_count = 10'000;

  std::string log_output_file;
  std::string json_output_file;
  easylog::Severity log_severity = easylog::Severity::DEBUG;
  bool log_console = false;

  static auto parse_from_args(int argc, char** argv) noexcept -> rfl::Result<WIMParams>;
};

struct VertexSet {
  std::vector<vertex_id_t> vertex_list;
  DynamicBitset mask;

  explicit VertexSet(vertex_id_t n) : vertex_list(), mask(n) {}

  explicit VertexSet(vertex_id_t n, std::span<vertex_id_t> vertices)
      : vertex_list(vertices.begin(), vertices.end()), mask(n) {
    for (auto v : vertices) {
      BOOST_ASSERT_MSG(v >= 0 && v < n, "Vertex index out of range [0, n)");
      mask.set(v);
    }
  }

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
  using Allocator = std::pmr::polymorphic_allocator<vertex_id_t>;

  static inline auto* allocator = new std::pmr::monotonic_buffer_resource{};

  const InvAdjacencyList<WIMEdge>* inv_graph;
  std::discrete_distribution<vertex_id_t> center_distribution;
  // sketches[i] = { v1, v2... }, all the vertices in the i-th RR-sketch
  std::pmr::vector<std::pmr::vector<vertex_id_t>> sketches;
  // inv_sketches[v] = { i1, i2 ... }, all the RR-sketches that contain v
  std::pmr::vector<std::pmr::vector<size_t>> inv_sketches;

  explicit RRSketchSet(const InvAdjacencyList<WIMEdge>* inv_graph, std::span<const vertex_weight_t> vertex_weights)
      : inv_graph(inv_graph), center_distribution(vertex_weights.begin(), vertex_weights.end()), sketches(allocator),
        inv_sketches(allocator) {
    BOOST_ASSERT_MSG(vertex_weights.size() == graph::num_vertices(*inv_graph),
                     "Mismatch between # of vertices in the graph and # of vertices in the weight list.");
    auto n = graph::num_vertices(*inv_graph);
    inv_sketches.assign(n, std::pmr::vector<size_t>{allocator});
  }

  auto num_vertices() const -> vertex_id_t {
    return static_cast<vertex_id_t>(inv_sketches.size());
  }

  auto num_sketches() const -> size_t {
    return sketches.size();
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

auto wim_experiment(const AdjacencyListPair<WIMEdge>& graph, const WIMParams& params) noexcept -> rfl::Result<json>;
