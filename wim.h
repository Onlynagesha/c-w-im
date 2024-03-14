#pragma once

#include "graph_types.h"
#include "utils/boost_assert.h"
#include "utils/dynamic_bitset.h"
#include "utils/graph.h"
#include <rfl/Size.hpp>
#include <rfl/Validator.hpp>
#include <rfl/comparisons.hpp>

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
      : VertexSet(n, std::span{vertices.begin(), vertices.end()}) {}

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

  auto n_vertices() const -> vertex_id_t {
    return static_cast<vertex_id_t>(inv_sketches.size());
  }

  auto n_sketches() const -> size_t {
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

auto wim_simulate(const AdjacencyList<WIMEdge>& graph, const VertexSet& seeds, uint64_t try_count) noexcept
    -> rfl::Result<double>;

auto wim_simulate_w(const AdjacencyList<WIMEdge>& graph, std::span<const vertex_weight_t> vertex_weights,
                    const VertexSet& seeds, uint64_t try_count) noexcept -> rfl::Result<double>;

auto wbim_simulate(const AdjacencyList<WBIMEdge>& graph, const VertexSet& seeds, const VertexSet& boosted_vertices,
                   uint64_t try_count) noexcept -> rfl::Result<double>;

auto wbim_simulate_w(const AdjacencyList<WBIMEdge>& graph, std::span<const vertex_weight_t> vertex_weights,
                     const VertexSet& seeds, const VertexSet& boosted_vertices, uint64_t try_count) noexcept
    -> rfl::Result<double>;

// result[v] = A value in range [0, 1]. Frequency that v received message in all the simulations
auto wim_simulate_p(const AdjacencyList<WIMEdge>& graph, const VertexSet& seeds, uint64_t try_count)
    -> rfl::Result<std::vector<double>>;

auto wbim_simulate_p(const AdjacencyList<WBIMEdge>& graph, const VertexSet& seeds, const VertexSet& boosted_vertices,
                     uint64_t try_count) -> rfl::Result<std::vector<double>>;

namespace details::wim {
template <class T>
auto as_vertex_set(vertex_id_t n, T&& vertices) -> decltype(auto) {
  if constexpr (std::is_same_v<std::remove_cvref_t<T>, VertexSet>) {
    return std::forward<T>(vertices);
  } else if constexpr (contiguous_range_of_exactly<std::remove_cvref_t<T>, vertex_id_t>) {
    return VertexSet(n, vertices);
  } else if constexpr (same_as_either<std::remove_cvref_t<T>, std::nullptr_t, std::nullopt_t>) {
    return VertexSet(n, {});
  } else {
    static_assert(rfl::always_false_v<T>, "Invalid type: Expects VertexSet or some contiguous range of vertex_id_t.");
  }
}
} // namespace details::wim

template <is_edge_property E, class VertexWeights, class SeedsOrList, class BoostedOrList>
inline auto simulate(const AdjacencyList<E>& graph, VertexWeights&& vertex_weights, SeedsOrList&& seeds_or_list,
                     BoostedOrList&& boosted_or_list, uint64_t try_count) -> rfl::Result<double> {
  constexpr auto USES_VERTEX_WEIGHTS =
      same_as_either<std::remove_cvref_t<VertexWeights>, std::nullptr_t, std::nullopt_t>;
  static_assert(USES_VERTEX_WEIGHTS || std::is_convertible_v<VertexWeights, std::span<vertex_id_t>>,
                "Invalid vertex weight type: Expects a contiguous range of vertex_id_t, "
                "or either of nullptr or std::nullopt as placeholder.");

  auto n = graph::num_vertices(graph);
  auto&& seeds = details::wim::as_vertex_set(n, seeds_or_list);
  if constexpr (std::is_same_v<E, WIMEdge>) {
    if (USES_VERTEX_WEIGHTS) {
      return wim_simulate_w(graph, vertex_weights, seeds, try_count);
    } else {
      return wim_simulate(graph, seeds, try_count);
    }
  } else if constexpr (std::is_same_v<E, WBIMEdge>) {
    auto&& boosted_vertices = details::wim::as_vertex_set(n, boosted_or_list);
    if (USES_VERTEX_WEIGHTS) {
      return wbim_simulate_w(graph, vertex_weights, seeds, boosted_vertices, try_count);
    } else {
      return wbim_simulate(graph, seeds, boosted_vertices, try_count);
    }
  } else {
    static_assert(rfl::always_false_v<E>, "Invalid edge type.");
  }
}

template <is_edge_property E, class SeedsOrList, class BoostedOrList>
inline auto simulate_p(const AdjacencyList<E>& graph, SeedsOrList&& seeds_or_list, BoostedOrList&& boosted_or_list,
                       uint64_t try_count) -> rfl::Result<std::vector<double>> {
  auto n = graph::num_vertices(graph);
  auto&& seeds = details::wim::as_vertex_set(n, seeds_or_list);
  if constexpr (std::is_same_v<E, WIMEdge>) {
    return wim_simulate_p(graph, seeds, try_count);
  } else if constexpr (std::is_same_v<E, WBIMEdge>) {
    auto&& boosted_vertices = details::wim::as_vertex_set(n, boosted_or_list);
    return wbim_simulate_p(graph, seeds, boosted_vertices, try_count);
  } else {
    static_assert(rfl::always_false_v<E>, "Invalid edge type.");
  }
}

struct DetectProbabilityFromSeedsResult {
  std::vector<edge_probability_t> p_in;
  std::vector<edge_probability_t> p_in_boosted;

  DetectProbabilityFromSeedsResult() = default;
  explicit DetectProbabilityFromSeedsResult(vertex_id_t n) : p_in(n, 0.0_ep), p_in_boosted(n, 0.0_ep) {}

  auto swap(DetectProbabilityFromSeedsResult& rhs) -> void {
    p_in.swap(rhs.p_in);
    p_in_boosted.swap(rhs.p_in_boosted);
  }

  auto equals_with(const DetectProbabilityFromSeedsResult& rhs) const -> bool {
    return ranges::equal(p_in, rhs.p_in) && ranges::equal(p_in_boosted, rhs.p_in_boosted);
  }
};

auto wim_detect_probability_from_seeds(const InvAdjacencyList<WIMEdge>& inv_graph, const VertexSet& seeds,
                                       vertex_id_t max_distance) -> rfl::Result<DetectProbabilityFromSeedsResult>;

auto wbim_detect_probability_from_seeds(const InvAdjacencyList<WBIMEdge>& inv_graph, const VertexSet& seeds,
                                        vertex_id_t max_distance) -> rfl::Result<DetectProbabilityFromSeedsResult>;

template <is_edge_property E>
auto detect_probability_from_seeds(const InvAdjacencyList<E>& inv_graph, const VertexSet& seeds,
                                   vertex_id_t max_distance) -> rfl::Result<DetectProbabilityFromSeedsResult> {
  if constexpr (std::is_same_v<E, WIMEdge>) {
    return wim_detect_probability_from_seeds(inv_graph, seeds, max_distance);
  } else if constexpr (std::is_same_v<E, WBIMEdge>) {
    return wbim_detect_probability_from_seeds(inv_graph, seeds, max_distance);
  } else {
    static_assert(rfl::always_false_v<E>, "Invalid edge type.");
  }
}
