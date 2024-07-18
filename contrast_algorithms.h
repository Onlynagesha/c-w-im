#pragma once

#include "graph_types.h"
#include "utils/result.h"
#include <rfl/Validator.hpp>
#include <rfl/comparisons.hpp>

auto wim_max_out_degree(const AdjacencyList<WIMEdge>& graph, vertex_id_t top_k = 0) -> std::vector<vertex_id_t>;

auto wbim_max_out_degree(const AdjacencyList<WBIMEdge>& graph, vertex_id_t top_k = 0) -> std::vector<vertex_id_t>;

auto wim_max_out_degree_i(const InvAdjacencyList<WIMEdge>& inv_graph, vertex_id_t top_k = 0)
    -> std::vector<vertex_id_t>;

auto wbim_max_out_degree_i(const InvAdjacencyList<WBIMEdge>& inv_graph, vertex_id_t top_k = 0)
    -> std::vector<vertex_id_t>;

template <int IsInv, is_edge_property E>
inline auto max_out_degree(const graph::adjacency<IsInv, E>& graph, vertex_id_t top_k = 0) {
  if constexpr (std::is_same_v<E, WIMEdge>) {
    if constexpr (IsInv) {
      return wim_max_out_degree_i(graph, top_k);
    } else {
      return wim_max_out_degree(graph, top_k);
    }
  } else if constexpr (std::is_same_v<E, WBIMEdge>) {
    if constexpr (IsInv) {
      return wbim_max_out_degree_i(graph, top_k);
    } else {
      return wbim_max_out_degree(graph, top_k);
    }
  } else {
    static_assert(!rfl::always_false_v<E>, "Invalid edge type.");
  }
}

auto wim_max_out_strength(const AdjacencyList<WIMEdge>& graph, vertex_id_t top_k = 0) -> std::vector<vertex_id_t>;

auto wbim_max_out_strength(const AdjacencyList<WBIMEdge>& graph, vertex_id_t top_k = 0) -> std::vector<vertex_id_t>;

auto wim_max_out_strength_i(const InvAdjacencyList<WIMEdge>& inv_graph, vertex_id_t top_k = 0)
    -> std::vector<vertex_id_t>;

auto wbim_max_out_strength_i(const InvAdjacencyList<WBIMEdge>& inv_graph, vertex_id_t top_k = 0)
    -> std::vector<vertex_id_t>;

template <int IsInv, is_edge_property E>
inline auto max_out_strength(const graph::adjacency<IsInv, E>& graph, vertex_id_t top_k = 0) {
  if constexpr (std::is_same_v<E, WIMEdge>) {
    if constexpr (IsInv) {
      return wim_max_out_strength_i(graph, top_k);
    } else {
      return wim_max_out_strength(graph, top_k);
    }
  } else if constexpr (std::is_same_v<E, WBIMEdge>) {
    if constexpr (IsInv) {
      return wbim_max_out_strength_i(graph, top_k);
    } else {
      return wbim_max_out_strength(graph, top_k);
    }
  } else {
    static_assert(!rfl::always_false_v<E>, "Invalid edge type.");
  }
}

struct PagerankParams {
  using DampingFactor = rfl::Validator<edge_probability_t, rfl::Minimum<0>, rfl::Maximum<1>>;
  DampingFactor damping_factor = 0.85_ep;

  using Epsilon = rfl::Validator<edge_probability_t, rfl::Minimum<0>>;
  Epsilon epsilon = 1e-6_ep;

  vertex_id_t k = 0;
  uint64_t n_iterations = 1'000;
  bool transpose = false;
  bool uses_vertex_weight = true;
  bool uses_edge_weight = true;
};

// vertex_weights = {} for unweighted vertices (which requires uses_vertex_weight == false)
auto wim_pagerank(const AdjacencyList<WIMEdge>& graph, const InvAdjacencyList<WIMEdge>& inv_graph,
                  std::span<const edge_probability_t> vertex_weights, const PagerankParams& params)
    -> rfl::Result<std::vector<edge_probability_t>>;

// vertex_weights = {} for unweighted vertices (which requires uses_vertex_weight == false)
auto wbim_pagerank(const AdjacencyList<WBIMEdge>& graph, const InvAdjacencyList<WBIMEdge>& inv_graph,
                   std::span<const edge_probability_t> vertex_weights, const PagerankParams& params)
    -> rfl::Result<std::vector<edge_probability_t>>;

template <is_edge_property E>
inline auto pagerank(const AdjacencyList<E>& graph, const InvAdjacencyList<E>& inv_graph,
                     std::span<const edge_probability_t> vertex_weights, const PagerankParams& params) {
  if constexpr (std::is_same_v<E, WIMEdge>) {
    return wim_pagerank(graph, inv_graph, vertex_weights, params);
  } else if constexpr (std::is_same_v<E, WBIMEdge>) {
    return wbim_pagerank(graph, inv_graph, vertex_weights, params);
  } else {
    static_assert(rfl::always_false_v<E>, "Invalid edge type.");
  }
}

template <is_edge_property E>
inline auto pagerank(const AdjacencyListPair<E>& graph, const PagerankParams& params) {
  const auto& [adj_list, inv_adj_list, vertex_weights] = graph;
  return pagerank(adj_list, inv_adj_list, vertex_weights, params);
}

auto wim_min_pagerank(const AdjacencyList<WIMEdge>& graph, const InvAdjacencyList<WIMEdge>& inv_graph,
                      std::span<const edge_probability_t> vertex_weights, const PagerankParams& params)
    -> rfl::Result<std::vector<vertex_id_t>>;

auto wim_max_pagerank(const AdjacencyList<WIMEdge>& graph, const InvAdjacencyList<WIMEdge>& inv_graph,
                      std::span<const edge_probability_t> vertex_weights, const PagerankParams& params)
    -> rfl::Result<std::vector<vertex_id_t>>;

auto wbim_min_pagerank(const AdjacencyList<WBIMEdge>& graph, const InvAdjacencyList<WBIMEdge>& inv_graph,
                       std::span<const edge_probability_t> vertex_weights, const PagerankParams& params)
    -> rfl::Result<std::vector<vertex_id_t>>;

auto wbim_max_pagerank(const AdjacencyList<WBIMEdge>& graph, const InvAdjacencyList<WBIMEdge>& inv_graph,
                       std::span<const edge_probability_t> vertex_weights, const PagerankParams& params)
    -> rfl::Result<std::vector<vertex_id_t>>;

template <same_as_either<ranges::less, ranges::greater> Compare, is_edge_property E>
inline auto best_pagerank(const AdjacencyList<E>& graph, const InvAdjacencyList<E>& inv_graph,
                          std::span<const edge_probability_t> vertex_weights, const PagerankParams& params)
    -> rfl::Result<std::vector<vertex_id_t>> {
  if constexpr (std::is_same_v<E, WIMEdge>) {
    if constexpr (std::is_same_v<Compare, ranges::less>) {
      return wim_min_pagerank(graph, inv_graph, vertex_weights, params);
    } else if constexpr (std::is_same_v<Compare, ranges::greater>) {
      return wim_max_pagerank(graph, inv_graph, vertex_weights, params);
    } else {
      static_assert(rfl::always_false_v<Compare>, "Invalid comparison policy.");
    }
  } else if constexpr (std::is_same_v<E, WBIMEdge>) {
    if constexpr (std::is_same_v<Compare, ranges::less>) {
      return wbim_min_pagerank(graph, inv_graph, vertex_weights, params);
    } else if constexpr (std::is_same_v<Compare, ranges::greater>) {
      return wbim_max_pagerank(graph, inv_graph, vertex_weights, params);
    } else {
      static_assert(rfl::always_false_v<Compare>, "Invalid comparison policy.");
    }
  } else {
    static_assert(rfl::always_false_v<E>, "Invalid edge type.");
  }
}

template <same_as_either<ranges::less, ranges::greater> Compare, is_edge_property E>
inline auto best_pagerank(const AdjacencyListPair<E>& graph, const PagerankParams& params) {
  const auto& [adj_list, inv_adj_list, vertex_weights] = graph;
  return best_pagerank<Compare>(adj_list, inv_adj_list, vertex_weights, params);
}

template <is_edge_property E>
inline auto max_pagerank(const AdjacencyList<E>& graph, const InvAdjacencyList<E>& inv_graph,
                         std::span<const edge_probability_t> vertex_weights, const PagerankParams& params)
    -> rfl::Result<std::vector<vertex_id_t>> {
  return best_pagerank<ranges::greater>(graph, inv_graph, vertex_weights, params);
}

template <is_edge_property E>
inline auto max_pagerank(const AdjacencyListPair<E>& graph, const PagerankParams& params) {
  const auto& [adj_list, inv_adj_list, vertex_weights] = graph;
  return max_pagerank(adj_list, inv_adj_list, vertex_weights, params);
}

template <is_edge_property E>
inline auto min_pagerank(const AdjacencyList<E>& graph, const InvAdjacencyList<E>& inv_graph,
                         std::span<const edge_probability_t> vertex_weights, const PagerankParams& params)
    -> rfl::Result<std::vector<vertex_id_t>> {
  return best_pagerank<ranges::less>(graph, inv_graph, vertex_weights, params);
}

template <is_edge_property E>
inline auto min_pagerank(const AdjacencyListPair<E>& graph, const PagerankParams& params) {
  const auto& [adj_list, inv_adj_list, vertex_weights] = graph;
  return min_pagerank(adj_list, inv_adj_list, vertex_weights, params);
}

struct IMRankParams {
  vertex_id_t k = 0;
  uint64_t n_iterations = 100;
  rfl::Validator<uint64_t, rfl::Minimum<1>> n_iterations_before_topk_fixed = 3;
};

// Returns top-k vertex indices
auto wim_imrank(const InvAdjacencyList<WIMEdge>& inv_graph, const IMRankParams& params) -> std::vector<vertex_id_t>;
auto wbim_imrank(const InvAdjacencyList<WBIMEdge>& inv_graph, const IMRankParams& params) -> std::vector<vertex_id_t>;

template <is_edge_property E>
inline auto imrank(const InvAdjacencyList<E>& inv_graph, const IMRankParams& params) {
  if constexpr (std::is_same_v<E, WIMEdge>) {
    return wim_imrank(inv_graph, params);
  } else if constexpr (std::is_same_v<E, WBIMEdge>) {
    return wbim_imrank(inv_graph, params);
  } else {
    static_assert(rfl::always_false_v<E>, "Invalid edge type.");
  }
}

struct RRobustSCCParams {
  uint64_t r;
};

auto wim_r_robust_scc(const AdjacencyList<WIMEdge>& graph, const InvAdjacencyList<WIMEdge>& inv_graph,
                      const RRobustSCCParams& params) -> rfl::Result<std::vector<vertex_id_t>>;