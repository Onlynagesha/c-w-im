#pragma once

#include "utils/graph.h"
#include "utils/result.h"
#include "utils/utils.h"
#include <compare>
#include <concepts>
#include <rfl/Result.hpp>

using vertex_weight_t = float;
using edge_probability_t = float;

inline constexpr auto operator""_vw(long double c) {
  return static_cast<vertex_weight_t>(c);
}

inline constexpr auto operator""_ep(long double w) {
  return static_cast<edge_probability_t>(w);
}

struct WIMEdge {
  float p;
  float p_seed;

  auto is_valid() const -> bool {
    return 0.0f <= p && p <= p_seed && p_seed <= 1.0f;
  }
  auto rand_test(bool is_seed) const -> bool {
    return rand_bool(is_seed ? p_seed : p);
  }
  constexpr auto operator<=>(const WIMEdge& rhs) const = default;
};

struct WBIMEdge {
  float p;
  float p_seed;
  float p_boost;

  auto is_valid() const -> bool {
    return 0.0f <= p && p <= std::min(p_seed, p_boost) <= std::max(p_seed, p_boost) <= 1.0f;
  }
  auto rand_test(bool is_seed, bool is_boosted) const -> bool {
    return rand_bool(is_seed ? p_seed : is_boosted ? p_boost : p);
  }
  constexpr auto operator<=>(const WBIMEdge& rhs) const = default;
};

template <class E>
concept is_edge_property = std::same_as<E, WIMEdge> || std::same_as<E, WBIMEdge>;

template <is_edge_property E>
struct ReadGraphResult {
  DirectedEdgeList<E> edge_list;
  std::vector<vertex_weight_t> vertex_weights;
};

using WIMReadGraphResult = ReadGraphResult<WIMEdge>;
using WBIMReadGraphResult = ReadGraphResult<WBIMEdge>;

auto read_directed_wim_edge_list(const std::string& input_file) noexcept -> rfl::Result<WIMReadGraphResult>;

auto read_directed_wbim_edge_list(const std::string& input_file) noexcept -> rfl::Result<WBIMReadGraphResult>;

auto write_directed_wim_edge_list(const DirectedEdgeList<WIMEdge>& graph,
                                  std::span<const vertex_weight_t> vertex_weights,
                                  const std::string& output_file) noexcept -> ResultVoid;

auto write_directed_wbim_edge_list(const DirectedEdgeList<WBIMEdge>& graph,
                                   std::span<const vertex_weight_t> vertex_weights,
                                   const std::string& output_file) noexcept -> ResultVoid;

inline auto write_directed_wim_edge_list_r(const WIMReadGraphResult& graph, const std::string& output_file) {
  return write_directed_wim_edge_list(graph.edge_list, graph.vertex_weights, output_file);
}

inline auto write_directed_wbim_edge_list_r(const WBIMReadGraphResult& graph, const std::string& output_file) {
  return write_directed_wbim_edge_list(graph.edge_list, graph.vertex_weights, output_file);
}
