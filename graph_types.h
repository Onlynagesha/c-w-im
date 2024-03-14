#pragma once

#include "utils/graph.h"
#include "utils/result.h"
#include "utils/utils.h"
#include <compare>
#include <concepts>
#include <rfl/Result.hpp>
#include <rfl/always_false.hpp>

using vertex_weight_t = float;
using edge_probability_t = float;

inline constexpr auto operator""_vw(long double c) {
  return static_cast<vertex_weight_t>(c);
}

inline constexpr auto operator""_ep(long double w) {
  return static_cast<edge_probability_t>(w);
}

struct WIMEdge {
  edge_probability_t p;
  edge_probability_t p_seed;

  auto is_valid() const -> bool {
    auto [p_min, p_max] = std::minmax(p, p_seed);
    return 0.0_ep <= p_min && p_max <= 1.0_ep;
  }
  auto rand_test(bool is_seed) const -> bool {
    return rand_bool(is_seed ? p_seed : p);
  }
  constexpr auto operator<=>(const WIMEdge& rhs) const = default;
};

struct WBIMEdge {
  edge_probability_t p;
  edge_probability_t p_boost;

  auto is_valid() const -> bool {
    auto [p_min, p_max] = std::minmax(p, p_boost);
    return 0.0_ep <= p_min && p_max <= 1.0_ep;
  }
  auto rand_test(bool is_boosted) const -> bool {
    return rand_bool(is_boosted ? p_boost : p);
  }
  constexpr auto operator<=>(const WBIMEdge& rhs) const = default;
};

template <class E>
concept is_edge_property = std::same_as<E, WIMEdge> || std::same_as<E, WBIMEdge>;

template <is_edge_property E>
inline auto get_p_seed_or_boost(const E& e) -> edge_probability_t {
  if constexpr (requires { e.p_seed; }) {
    return e.p_seed;
  } else if constexpr (requires { e.p_boost; }) {
    return e.p_boost;
  } else {
    return e.p;
  }
}

template <is_edge_property E>
inline auto get_p_seed(const E& e) -> edge_probability_t {
  if constexpr (requires { e.p_seed; }) {
    return e.p_seed;
  } else {
    return e.p;
  }
}

template <is_edge_property E>
inline auto get_p_boost(const E& e) -> edge_probability_t {
  if constexpr (requires { e.p_boost; }) {
    return e.p_boost;
  } else {
    return e.p;
  }
};

template <is_edge_property E>
struct ReadGraphResult {
  DirectedEdgeList<E> edge_list;
  std::vector<vertex_weight_t> vertex_weights;
};

using WIMReadGraphResult = ReadGraphResult<WIMEdge>;
using WBIMReadGraphResult = ReadGraphResult<WBIMEdge>;

template <is_edge_property E>
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

  auto n_vertices() const {
    return graph::num_vertices(adj_list);
  }

  auto n_edges() const {
    return adj_list.num_edges();
  }

  auto graph_n_m() const {
    return std::tuple{n_vertices(), n_edges()};
  }
};

using WIMAdjacencyListPair = AdjacencyListPair<WIMEdge>;
using WBIMAdjacencyListPair = AdjacencyListPair<WBIMEdge>;

auto read_directed_wim_edge_list(const std::string& input_file) noexcept -> rfl::Result<WIMReadGraphResult>;

auto read_directed_wbim_edge_list(const std::string& input_file) noexcept -> rfl::Result<WBIMReadGraphResult>;

auto read_directed_wim_adjacency_lists(const std::string& input_file) noexcept -> rfl::Result<WIMAdjacencyListPair>;

auto read_directed_wbim_adjacency_lists(const std::string& input_file) noexcept -> rfl::Result<WBIMAdjacencyListPair>;

template <is_edge_property E>
inline auto read_directed_edge_list(const std::string& input_file) {
  if constexpr (std::is_same_v<E, WIMEdge>) {
    return read_directed_wim_edge_list(input_file);
  } else if constexpr (std::is_same_v<E, WBIMEdge>) {
    return read_directed_wbim_edge_list(input_file);
  } else {
    static_assert(rfl::always_false_v<E>, "Invalid edge type.");
  }
}

template <is_edge_property E>
inline auto read_directed_adjacency_lists(const std::string& input_file) {
  if constexpr (std::is_same_v<E, WIMEdge>) {
    return read_directed_wim_adjacency_lists(input_file);
  } else if constexpr (std::is_same_v<E, WBIMEdge>) {
    return read_directed_wbim_adjacency_lists(input_file);
  } else {
    static_assert(rfl::always_false_v<E>, "Invalid edge type.");
  }
}

auto write_directed_wim_edge_list(const DirectedEdgeList<WIMEdge>& graph,
                                  std::span<const vertex_weight_t> vertex_weights,
                                  const std::string& output_file) noexcept -> ResultVoid;

auto write_directed_wbim_edge_list(const DirectedEdgeList<WBIMEdge>& graph,
                                   std::span<const vertex_weight_t> vertex_weights,
                                   const std::string& output_file) noexcept -> ResultVoid;

template <is_edge_property E>
inline auto write_directed_edge_list(const DirectedEdgeList<E>& graph, std::span<const vertex_weight_t> vertex_weights,
                                     const std::string& output_file) {
  if constexpr (std::is_same_v<E, WIMEdge>) {
    return write_directed_wim_edge_list(graph, vertex_weights, output_file);
  } else if constexpr (std::is_same_v<E, WBIMEdge>) {
    return write_directed_wbim_edge_list(graph, vertex_weights, output_file);
  } else {
    static_assert(rfl::always_false_v<E>, "Invalid edge type.");
  }
}

template <is_edge_property E>
inline auto write_directed_edge_list_r(const ReadGraphResult<E>& graph, const std::string& output_file) {
  return write_directed_edge_list(graph.edge_list, graph.vertex_weights, output_file);
}
