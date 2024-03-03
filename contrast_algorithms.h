#pragma once

#include "graph_types.h"
#include "utils/result.h"
#include <rfl/Validator.hpp>
#include <rfl/comparisons.hpp>

struct PagerankParams {
  using DampingFactor = rfl::Validator<edge_probability_t, rfl::Minimum<0>, rfl::Maximum<1>>;
  DampingFactor damping_factor = 0.85;
  uint64_t n_iterations = 1'000;
  bool uses_vertex_weight = true;
  bool uses_edge_weight = true;
};

auto wim_pagerank(const AdjacencyList<WIMEdge>& graph, const InvAdjacencyList<WIMEdge>& inv_graph,
                  std::span<const edge_probability_t> vertex_weights, const PagerankParams& params)
    -> rfl::Result<std::vector<edge_probability_t>>;
