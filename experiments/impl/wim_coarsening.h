#pragma once

#include "experiments.h"
#include "experiments/internal_states.h"

// Results are written to json_root
namespace exp_impl {
// Expands to the original graph -> Simulation
auto do_wim_experiment_expand_simulate(const AdjacencyListPair<WIMEdge>& original_graph,
                                       std::span<const WIMCoarsenGraphBriefResult> coarsen_results,
                                       exp_states::VertexListList coarsened_seed_lists, uint64_t simulation_try_count,
                                       const WIMSketchingParams& sketching_params,
                                       const MultiLevelParams& multi_level_params, json* json_root) -> ResultVoid;

// Gets seeds from the coarsened graph -> Expands to the original graph -> Simulation
auto do_wim_experiment_with_coarsened(const AdjacencyListPair<WIMEdge>& original_graph,
                                      std::span<const WIMCoarsenGraphBriefResult> coarsen_results,
                                      const CommonExperimentParams& common_params,
                                      const WIMSketchingParams& sketching_params,
                                      const MultiLevelParams& multi_level_params, json* json_root) -> ResultVoid;

// Overload with different parameter types.
inline auto do_wim_experiment_with_coarsened(const AdjacencyListPair<WIMEdge>& original_graph,
                                             std::span<const WIMCoarsenGraphBriefResult> coarsen_results,
                                             const WIMCoarseningExperimentParams& params, json* json_root)
    -> ResultVoid {
  return do_wim_experiment_with_coarsened( //
      original_graph, coarsen_results, *params.common, *params.sketching, *params.multi_level, json_root);
}

// Gets boosted vertices from the coarsened graph -> Expands to the original graph -> Simulation
auto do_wbim_experiment_with_coarsened(const AdjacencyListPair<WBIMEdge>& original_graph, const VertexSet& seeds,
                                       std::span<const WBIMCoarsenGraphBriefResult> coarsen_results,
                                       const CommonExperimentParams& common_params,
                                       const WBIMSketchingParams& sketching_params,
                                       const MultiLevelParams& multi_level_params, json* json_root) -> ResultVoid;

// Overload with different parameter types.
inline auto do_wbim_experiment_with_coarsened(const AdjacencyListPair<WBIMEdge>& original_graph, const VertexSet& seeds,
                                              std::span<const WBIMCoarsenGraphBriefResult> coarsen_results,
                                              const WBIMCoarseningExperimentParams& params, json* json_root)
    -> ResultVoid {
  return do_wbim_experiment_with_coarsened( //
      original_graph, seeds, coarsen_results, *params.common, *params.sketching, *params.multi_level, json_root);
}

auto do_wim_coarsening_experiment(const AdjacencyListPair<WIMEdge>& graph, const WIMCoarseningExperimentParams& params,
                                  json* json_root) -> ResultVoid;

auto do_wbim_coarsening_experiment(const AdjacencyListPair<WBIMEdge>& graph, const VertexSet& seeds,
                                   const WBIMCoarseningExperimentParams& params, json* json_root) -> ResultVoid;
} // namespace exp_impl
