#pragma once

#include "coarsening_details_types.h"
#include "graph_types.h"
#include "utils/result.h"
#include <array>
#include <rfl/Rename.hpp>
#include <rfl/Validator.hpp>
#include <rfl/comparisons.hpp>

enum class NeighborMatchRule { HEM_P_MAX, HEM_P_PRODUCT, LEM_P_MAX, LEM_P_PRODUCT };

enum class EdgeWeightRule { SEPARATE_SIMPLE, MERGED_SIMPLE, SEPARATE_PRECISE, MERGED_PRECISE };

enum class SeedEdgeWeightRule { AVERAGE, MAX, BEST_SEED_INDEX };

enum class BoostedEdgeWeightRule { AVERAGE, MAX, BEST_BOOSTED_INDEX };

enum class InOutHeuristicRule { UNIT, COUNT, P, W, SEED_P, SEED_W };

enum class BoostedSelectionRule { AS_SOURCE, AS_TARGET };

enum class VertexWeightRule { AVERAGE, AVERAGE_BY_PATHS, SUM };

enum class VertexExpandingRule { RANDOM, LOCAL, SIMULATIVE, ITERATIVE };

inline constexpr auto rule_prefix_is_seed(InOutHeuristicRule rule) {
  return rule == InOutHeuristicRule::SEED_P || rule == InOutHeuristicRule::SEED_W;
}

struct CoarseningParams {
  NeighborMatchRule neighbor_match_rule = NeighborMatchRule::HEM_P_MAX;
  EdgeWeightRule edge_weight_rule = EdgeWeightRule::SEPARATE_PRECISE;
  SeedEdgeWeightRule seed_edge_weight_rule = SeedEdgeWeightRule::BEST_SEED_INDEX;
  BoostedEdgeWeightRule boosted_edge_weight_rule = BoostedEdgeWeightRule::BEST_BOOSTED_INDEX;
  BoostedSelectionRule boosted_selection_rule = BoostedSelectionRule::AS_TARGET;
  InOutHeuristicRule in_out_heuristic_rule = InOutHeuristicRule::P;
  VertexWeightRule vertex_weight_rule = VertexWeightRule::AVERAGE_BY_PATHS;
  vertex_id_t max_distance_from_seed = 6;
};

struct ExpandingParams {
  VertexExpandingRule vertex_expanding_rule;
  // Used for MERGED only
  using NIterations = rfl::Rename<              //
      "expanding_n_iterations",                 //
      rfl::Validator<uint64_t, rfl::Minimum<1>> //
      >;
  NIterations n_iterations = 10;

  using SimulationTryCount = rfl::Rename<       //
      "expanding_simulation_try_count",         //
      rfl::Validator<uint64_t, rfl::Minimum<1>> //
      >;
  SimulationTryCount simulation_try_count = 1'000;
};

// ---- Step 1: Converts to undirected edge ----

auto merge_wim_edge_to_undirected(const AdjacencyList<WIMEdge>& graph, const CoarseningParams& params) noexcept
    -> AdjacencyList<edge_probability_t>;

auto merge_wbim_edge_to_undirected(const AdjacencyList<WBIMEdge>& graph, const CoarseningParams& params) noexcept
    -> AdjacencyList<edge_probability_t>;

template <is_edge_property E>
inline auto merge_edge_to_undirected(const AdjacencyList<E>& graph, const CoarseningParams& params) noexcept {
  if constexpr (std::is_same_v<E, WIMEdge>) {
    return merge_wim_edge_to_undirected(graph, params);
  } else if constexpr (std::is_same_v<E, WBIMEdge>) {
    return merge_wbim_edge_to_undirected(graph, params);
  } else {
    static_assert(rfl::always_false_v<E>, "Invalid edge type.");
  }
}

// ---- Step 2: Performs Mongoose matching algorithm ----

struct MongooseMatchResult {
  vertex_id_t n_groups;
  std::vector<vertex_id_t> group_id;
};

auto mongoose_match(const AdjacencyList<edge_probability_t>& bidir_graph, const CoarseningParams& params) noexcept
    -> MongooseMatchResult;

auto mongoose_match_s(const AdjacencyList<edge_probability_t>& bidir_graph, std::span<const vertex_id_t> seeds,
                      const CoarseningParams& params) noexcept -> MongooseMatchResult;

// ---- Step 3: Performs coarsening by the groups obtained previously ----

auto coarsen_wim_graph_by_match(const AdjacencyList<WIMEdge>& graph, const InvAdjacencyList<WIMEdge>& inv_graph,
                                std::span<const vertex_weight_t> vertex_weights, vertex_id_t n_groups,
                                std::span<const vertex_id_t> group_id, const CoarseningParams& params) noexcept
    -> rfl::Result<WIMCoarsenGraphBriefResult>;

auto coarsen_wim_graph_by_match_d(const AdjacencyList<WIMEdge>& graph, const InvAdjacencyList<WIMEdge>& inv_graph,
                                  std::span<const vertex_weight_t> vertex_weights, vertex_id_t n_groups,
                                  std::span<const vertex_id_t> group_id, const CoarseningParams& params) noexcept
    -> rfl::Result<WIMCoarsenGraphDetailedResult>;

auto coarsen_wbim_graph_by_match(const AdjacencyList<WBIMEdge>& graph, const InvAdjacencyList<WBIMEdge>& inv_graph,
                                 std::span<const vertex_weight_t> vertex_weights, const VertexSet& seeds,
                                 vertex_id_t n_groups, std::span<const vertex_id_t> group_id,
                                 const CoarseningParams& params) noexcept -> rfl::Result<WBIMCoarsenGraphBriefResult>;

auto coarsen_wbim_graph_by_match_d(const AdjacencyList<WBIMEdge>& graph, const InvAdjacencyList<WBIMEdge>& inv_graph,
                                   std::span<const vertex_weight_t> vertex_weights, const VertexSet& seeds,
                                   vertex_id_t n_groups, std::span<const vertex_id_t> group_id,
                                   const CoarseningParams& params) noexcept
    -> rfl::Result<WBIMCoarsenGraphDetailedResult>;

// ---- Combination of Step 1 - 3 ----

auto coarsen_wim_graph(const AdjacencyList<WIMEdge>& graph, const InvAdjacencyList<WIMEdge>& inv_graph,
                       std::span<const vertex_weight_t> vertex_weights, const CoarseningParams& params) noexcept
    -> rfl::Result<WIMCoarsenGraphBriefResult>;

inline auto coarsen_wim_graph_p(const WIMAdjacencyListPair& graph, const CoarseningParams& params) noexcept
    -> rfl::Result<WIMCoarsenGraphBriefResult> {
  return coarsen_wim_graph(graph.adj_list, graph.inv_adj_list, graph.vertex_weights, params);
}

auto coarsen_wim_graph_d(const AdjacencyList<WIMEdge>& graph, const InvAdjacencyList<WIMEdge>& inv_graph,
                         std::span<const vertex_weight_t> vertex_weights, const CoarseningParams& params) noexcept
    -> rfl::Result<WIMCoarsenGraphDetailedResult>;

inline auto coarsen_wim_graph_d_p(const WIMAdjacencyListPair& graph, const CoarseningParams& params) noexcept
    -> rfl::Result<WIMCoarsenGraphDetailedResult> {
  return coarsen_wim_graph_d(graph.adj_list, graph.inv_adj_list, graph.vertex_weights, params);
}

auto coarsen_wbim_graph(const AdjacencyList<WBIMEdge>& graph, const InvAdjacencyList<WBIMEdge>& inv_graph,
                        std::span<const vertex_weight_t> vertex_weights, const VertexSet& seeds,
                        const CoarseningParams& params) noexcept -> rfl::Result<WBIMCoarsenGraphBriefResult>;

inline auto coarsen_wbim_graph_p(const WBIMAdjacencyListPair& graph, const VertexSet& seeds,
                                 const CoarseningParams& params) noexcept -> rfl::Result<WBIMCoarsenGraphBriefResult> {
  return coarsen_wbim_graph(graph.adj_list, graph.inv_adj_list, graph.vertex_weights, seeds, params);
}

auto coarsen_wbim_graph_d(const AdjacencyList<WBIMEdge>& graph, const InvAdjacencyList<WBIMEdge>& inv_graph,
                          std::span<const vertex_weight_t> vertex_weights, const VertexSet& seeds,
                          const CoarseningParams& params) noexcept -> rfl::Result<WBIMCoarsenGraphDetailedResult>;

inline auto coarsen_wbim_graph_d_p(const WBIMAdjacencyListPair& graph, const VertexSet& seeds,
                                   const CoarseningParams& params) noexcept
    -> rfl::Result<WBIMCoarsenGraphDetailedResult> {
  return coarsen_wbim_graph_d(graph.adj_list, graph.inv_adj_list, graph.vertex_weights, seeds, params);
}

// ---- Step 4: Expanding seeds ----

struct SelectBestSeedResult {
  size_t index_in_group;
};

struct SelectBestBoostedResult {
  size_t index_in_group;
};

auto select_best_seed_in_group(const WIMCoarsenedVertexDetails& v) noexcept -> SelectBestSeedResult;

auto select_best_boosted_in_group(const WBIMCoarsenedVertexDetails& v, BoostedSelectionRule rule) noexcept
    -> SelectBestBoostedResult;

struct ExpandSeedResult {
  std::vector<vertex_id_t> expanded_seeds;
};

struct ExpandBoostedResult {
  std::vector<vertex_id_t> expanded_boosted;
};

auto expand_wim_seed_vertices(const AdjacencyList<WIMEdge>& graph, std::span<const vertex_weight_t> vertex_weights,
                              const WIMCoarsenGraphBriefResult& coarsening_result,
                              std::span<const vertex_id_t> coarsened_seeds, const ExpandingParams& params) noexcept
    -> rfl::Result<ExpandSeedResult>;

inline auto further_expand_wim_seed_vertices(const WIMCoarsenGraphBriefResult& last_result,
                                             const WIMCoarsenGraphBriefResult& cur_result,
                                             std::span<const vertex_id_t> coarsened_seeds,
                                             const ExpandingParams& params) -> rfl::Result<ExpandSeedResult> {
  return expand_wim_seed_vertices(last_result.coarsened.adj_list, last_result.coarsened.vertex_weights, //
                                  cur_result, coarsened_seeds, params);
}

auto expand_wim_seed_vertices_d(const AdjacencyList<WIMEdge>& graph, std::span<const vertex_weight_t> vertex_weights,
                                const WIMCoarsenGraphDetailedResult& coarsening_result,
                                std::span<const vertex_id_t> coarsened_seeds, const ExpandingParams& params) noexcept
    -> rfl::Result<ExpandSeedResult>;

inline auto further_expand_wim_seed_vertices_d(const WIMCoarsenGraphDetailedResult& last_result,
                                               const WIMCoarsenGraphDetailedResult& cur_result,
                                               std::span<const vertex_id_t> coarsened_seeds,
                                               const ExpandingParams& params) -> rfl::Result<ExpandSeedResult> {
  return expand_wim_seed_vertices_d(last_result.coarsened.adj_list, last_result.coarsened.vertex_weights, //
                                    cur_result, coarsened_seeds, params);
}

auto expand_wbim_boosted_vertices(const AdjacencyList<WBIMEdge>& graph, std::span<const vertex_weight_t> vertex_weights,
                                  const WBIMCoarsenGraphBriefResult& coarsening_result,
                                  std::span<const vertex_id_t> coarsened_boosted,
                                  const ExpandingParams& params) noexcept -> rfl::Result<ExpandBoostedResult>;

inline auto further_expand_wbim_boosted_vertices(const WBIMCoarsenGraphBriefResult& last_result,
                                                 const WBIMCoarsenGraphBriefResult& cur_result,
                                                 std::span<const vertex_id_t> coarsened_boosted,
                                                 const ExpandingParams& params) -> rfl::Result<ExpandBoostedResult> {
  return expand_wbim_boosted_vertices(last_result.coarsened.adj_list, last_result.coarsened.vertex_weights, //
                                      cur_result, coarsened_boosted, params);
}

auto expand_wbim_boosted_vertices_d(const AdjacencyList<WBIMEdge>& graph,
                                    std::span<const vertex_weight_t> vertex_weights,
                                    const WBIMCoarsenGraphDetailedResult& coarsening_result,
                                    std::span<const vertex_id_t> coarsened_boosted,
                                    const ExpandingParams& params) noexcept -> rfl::Result<ExpandBoostedResult>;

inline auto further_expand_wbim_boosted_vertices_d(const WBIMCoarsenGraphDetailedResult& last_result,
                                                   const WBIMCoarsenGraphDetailedResult& cur_result,
                                                   std::span<const vertex_id_t> coarsened_boosted,
                                                   const ExpandingParams& params) -> rfl::Result<ExpandBoostedResult> {
  return expand_wbim_boosted_vertices_d(last_result.coarsened.adj_list, last_result.coarsened.vertex_weights, //
                                        cur_result, coarsened_boosted, params);
}
