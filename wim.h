#pragma once

#include "graph_types.h"
#include "utils/boost_assert.h"
#include "utils/flat_set.h"
#include "utils/graph.h"
#include "vertex_set.h"
#include <rfl/Size.hpp>
#include <rfl/Validator.hpp>
#include <rfl/comparisons.hpp>

#define WIM_DUMP_TYPES(F) \
  F(PRRSketch)            \
  F(PRRSketchSet)

#define DECLARE_FREE_DUMP_FUNCTION_FOR_WIM_TYPES(Type) \
  auto dump(const Type& obj, int indent = 0, int level = 0) noexcept -> std::string;

enum class PRRSketchEdgeState { BLOCKED, LIVE, LIVE_UPON_BOOST };

inline auto get_random_edge_state(const WBIMEdge& e) -> PRRSketchEdgeState {
  BOOST_ASSERT(e.is_valid());
  auto r = rand_float();
  if (r < e.p) {
    return PRRSketchEdgeState::LIVE;
  } else if (r < e.p_boost) {
    return PRRSketchEdgeState::LIVE_UPON_BOOST;
  } else {
    return PRRSketchEdgeState::BLOCKED;
  }
}

struct RRSketchSet {
  const InvAdjacencyList<WIMEdge>* inv_graph;
  std::discrete_distribution<vertex_id_t> center_distribution;
  // sketches[i] = { v1, v2... }, all the vertices in the i-th RR-sketch
  std::vector<std::vector<vertex_id_t>> sketches;
  // inv_sketches[v] = { i1, i2 ... }, all the RR-sketches that contain v
  std::vector<std::vector<size_t>> inv_sketches;

  RRSketchSet(const InvAdjacencyList<WIMEdge>* inv_graph, std::span<const vertex_weight_t> vertex_weights)
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

  auto rr_sketch_total_size_bytes() const -> size_t {
    return estimated_memory_usage(sketches) + estimated_memory_usage(inv_sketches);
  }

  auto rr_sketch_total_size_str() const -> std::string {
    return size_bytes_to_memory_str(rr_sketch_total_size_bytes());
  }

  // Appends r new RR-sketches
  auto append(size_t n_sketches) noexcept -> void;

  // Selects k seed vertices
  auto select(vertex_id_t k) const noexcept -> std::vector<vertex_id_t>;

  auto dump() const noexcept -> std::string;

private:
  // Appends a new RR-sketch
  auto append_single(std::span<vertex_id_t> vertices) noexcept -> void;
};

struct PRRSketch {
  static constexpr auto NULL_INDEX = std::numeric_limits<vertex_id_t>::max();
  static constexpr auto SUPER_SEED_INDEX = 0_vid;

  static_assert(to_signed(NULL_INDEX) < 0,
                "Assumption violated: as_signed(x) <= 0 is equivalent to x == 0 || x >= NULL_INDEX");

  // Seed vertices BEFORE MAPPING
  // (note: seeds are removed and then merged to a single "super seed" with index = 0)
  FlatSet<vertex_id_t> vertices;
  // Critical vertices (i.e. gain can be reached solely by setting the single vertex as boosted) BEFORE MAPPING
  FlatSet<vertex_id_t> critical_vertices;
  // The center vertex of current PRR-sketch, index BEFORE MAPPING.
  vertex_id_t center;
  // The center vertex of current PRR-sketch, index AFTER MAPPING.
  vertex_id_t mapped_center;
  // The PRR-sketch graph with edge states preserved, indices AFTER MAPPING, #0 as super seed.
  AdjacencyList<PRRSketchEdgeState> mapped_graph;
  // Transpose of mapped_graph, indices AFTER MAPPING
  InvAdjacencyList<PRRSketchEdgeState> inv_mapped_graph;

  auto mapped_vertex_index(vertex_id_t v) const -> vertex_id_t {
    auto it = vertices.find(v);
    if (it == vertices.end()) {
      return NULL_INDEX;
    }
    // 1 : Mapped indices start from 1 since 0 is reserved for the "super seed"
    return static_cast<vertex_id_t>(it - vertices.begin()) + 1;
  }

  auto dump(int indent = 0, int level = 0) const noexcept -> std::string;
};

struct PRRSketchSet {
  static constexpr auto NULL_INDEX = PRRSketch::NULL_INDEX;
  static constexpr auto SUPER_SEED_INDEX = PRRSketch::SUPER_SEED_INDEX;

  const AdjacencyList<WBIMEdge>* graph;
  const InvAdjacencyList<WBIMEdge>* inv_graph;
  const VertexSet* seeds;
  std::discrete_distribution<vertex_id_t> center_distribution;
  // sketches[i].critical_vertices = All the critical vertices in the i-th PRR-sketch
  std::vector<PRRSketch> sketches;
  // For each i in inv_critical[v], v is a critical vertex in the i-th PRR-sketch
  std::vector<std::vector<size_t>> inv_critical;
  // # of sketching attempts (success + failure)
  uint64_t total_n_attempts;

  PRRSketchSet(const AdjacencyList<WBIMEdge>* graph, const InvAdjacencyList<WBIMEdge>* inv_graph,
               std::span<const vertex_weight_t> vertex_weights, const VertexSet* seeds)
      : graph(graph), inv_graph(inv_graph), seeds(seeds),
        center_distribution(vertex_weights.begin(), vertex_weights.end()), //
        sketches(), inv_critical(), total_n_attempts(0) {
    auto [n, m] = graph_n_m(*graph);
    auto [ni, mi] = graph_n_m(*inv_graph);
    BOOST_ASSERT_MSG(n == ni && m == mi, //
                     "Mismatch between the given graph and the transpose graph.");
    BOOST_ASSERT_MSG(vertex_weights.size() == n, //
                     "Mismatch between # of vertices in the graph and # of vertices in the weight list.");
    inv_critical.assign(n, std::vector<size_t>{});
  }

  PRRSketchSet(const AdjacencyListPair<WBIMEdge>* graph, const VertexSet* seeds)
      : PRRSketchSet(&graph->adj_list, &graph->inv_adj_list, graph->vertex_weights, seeds) {}

  auto n_vertices() const -> vertex_id_t {
    return graph::num_vertices(*graph);
  }

  auto n_sketches() const -> size_t {
    return sketches.size();
  }

  auto sketch_n_vertices() const {
    return sketches | TRANSFORM_VIEW(graph::num_vertices(_1.mapped_graph));
  }

  auto sketch_n_edges() const {
    return sketches | TRANSFORM_VIEW(_1.mapped_graph.num_edges());
  }

  auto sketch_size_bytes() const {
    constexpr auto to_estimated_size_bytes = [&](const PRRSketch& sketch) -> size_t {
      auto res = sizeof(sketch.center) + sizeof(sketch.mapped_center);
      res += sketch.vertices.capacity() * sizeof(vertex_id_t);
      res += sketch.critical_vertices.capacity() * sizeof(vertex_id_t);
      // 2 : mapped + inv-mapped
      auto [n, m] = graph_n_m(sketch.mapped_graph);
      res += 2 * (sizeof(vertex_id_t) * n + (sizeof(vertex_id_t) + sizeof(PRRSketchEdgeState)) * m);
      return res;
    };
    return sketches | views::transform(to_estimated_size_bytes);
  }

  auto sketch_total_size_bytes() const -> size_t {
    return accumulate_sum(sketch_size_bytes());
  }

  auto average_sketch_n_vertices() const {
    BOOST_ASSERT_MSG(!sketches.empty(), "Requires at least 1 PRR-sketch to exist on calculating average size.");
    return 1.0 * accumulate_sum(sketch_n_vertices()) / sketches.size();
  }

  auto average_sketch_n_edges() const {
    BOOST_ASSERT_MSG(!sketches.empty(), "Requires at least 1 PRR-sketch to exist on calculating average size.");
    return 1.0 * accumulate_sum(sketch_n_edges()) / sketches.size();
  }

  auto average_sketch_size_bytes() const -> double {
    BOOST_ASSERT_MSG(!sketches.empty(), "Requires at least 1 PRR-sketch to exist on calculating average size.");
    return 1.0 * sketch_total_size_bytes() / sketches.size();
  }

  auto success_rate() const -> double {
    BOOST_ASSERT_MSG(total_n_attempts > 0, "Requires at least 1 attempts.");
    return 1.0 * sketches.size() / total_n_attempts;
  }

  // Appends r new PRR-sketches
  auto append(size_t n_sketches, vertex_id_t k) noexcept -> void;

  // Selects k boosted vertices
  auto select(vertex_id_t k) const noexcept -> std::vector<vertex_id_t>;

  // Selects k boosted vertices by critical vertices only (which is much faster)
  auto select_by_critical(vertex_id_t k) const noexcept -> std::vector<vertex_id_t>;

  auto dump(int indent = 0, int level = 0) const noexcept -> std::string;

private:
  // Appends a new PRR-sketch of given center.
  auto append_single(vertex_id_t center, vertex_id_t k, void* cache_raw_ptr) noexcept -> bool;
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
  constexpr auto NO_VERTEX_WEIGHTS = same_as_either<std::remove_cvref_t<VertexWeights>, std::nullptr_t, std::nullopt_t>;
  static_assert(NO_VERTEX_WEIGHTS || std::is_convertible_v<VertexWeights, std::span<const vertex_weight_t>>,
                "Invalid vertex weight type: Expects a contiguous range of vertex_id_t, "
                "or either of nullptr or std::nullopt as placeholder.");

  auto n = graph::num_vertices(graph);
  auto&& seeds = details::wim::as_vertex_set(n, seeds_or_list);
  if constexpr (std::is_same_v<E, WIMEdge>) {
    if (!NO_VERTEX_WEIGHTS) {
      return wim_simulate_w(graph, vertex_weights, seeds, try_count);
    } else {
      return wim_simulate(graph, seeds, try_count);
    }
  } else if constexpr (std::is_same_v<E, WBIMEdge>) {
    auto&& boosted_vertices = details::wim::as_vertex_set(n, boosted_or_list);
    if (!NO_VERTEX_WEIGHTS) {
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

struct WBIMActivationProbability {
  std::vector<edge_probability_t> p_in;
  std::vector<edge_probability_t> p_in_boosted;

  WBIMActivationProbability() = default;
  explicit WBIMActivationProbability(vertex_id_t n) : p_in(n, 0.0_ep), p_in_boosted(n, 0.0_ep) {}

  auto swap(WBIMActivationProbability& rhs) -> void {
    p_in.swap(rhs.p_in);
    p_in_boosted.swap(rhs.p_in_boosted);
  }

  auto equals_with(const WBIMActivationProbability& rhs) const -> bool {
    return ranges::equal(p_in, rhs.p_in) && ranges::equal(p_in_boosted, rhs.p_in_boosted);
  }
};

auto wbim_activation_probability_from_seeds(const InvAdjacencyList<WBIMEdge>& inv_graph, const VertexSet& seeds,
                                            vertex_id_t max_distance) -> rfl::Result<WBIMActivationProbability>;

WIM_DUMP_TYPES(DECLARE_FREE_DUMP_FUNCTION_FOR_WIM_TYPES)

#undef DECLARE_FREE_DUMP_FUNCTION_FOR_WIM_TYPES
