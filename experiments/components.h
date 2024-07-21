#pragma once

#include "experiments.h"
#include "experiments/internal_states.h"

namespace exp_components {
using exp_states::VertexList;
using exp_states::VertexListList;

// ---- Group A: Seed Generation ----

// Generates seeds for WBIM experiments
auto generate_wbim_seeds(const WBIMAdjacencyListPair& graph, const WBIMSeedGeneratingParams& params)
    -> rfl::Result<VertexSet>;

struct GetSeedsParams {
  HistogramShape histogram_shape;
};

struct GetBoostedParams {
  HistogramShape histogram_shape;
};

// ---- Group B: Seed / Boosted Vertices Selection ----

// Multiple experiments combined.
// Let [R1, R2 ... Rn] = sketching_params.n_sketches, Km = max in sketching_params.n_seeds,
// the function gets Km seed lists, which are obtained from R1 RR-sketches, R2 RR-sketches, etc.
auto do_wim_experiment_get_seeds(const WIMAdjacencyListPair& graph, const WIMSketchingParams& sketching_params,
                                 const GetSeedsParams& other_params)
    -> rfl::Result<exp_states::WIMSketchingGetSeedsResult>;

// Overload of different parameter type.
inline auto do_wim_experiment_get_seeds(const WIMAdjacencyListPair& graph, const CommonExperimentParams& common,
                                        const WIMSketchingParams& sketching)
    -> rfl::Result<exp_states::WIMSketchingGetSeedsResult> {
  return exp_components::do_wim_experiment_get_seeds(graph, sketching, {.histogram_shape = common.histogram_shape()});
}

// Multiple experiments combined.
// Let [R1, R2 ... Rn] = sketching_params.n_sketches, Km = max in sketching_params.n_boosted,
// the function gets Km boosted lists, which are obtained from R1 PRR-sketches, R2 PRR-sketches, etc.
auto do_wbim_experiment_get_boosted(const WBIMAdjacencyListPair& graph, const VertexSet& seeds,
                                    const WBIMSketchingParams& sketching_params, const GetBoostedParams& other_params)
    -> rfl::Result<exp_states::WBIMSketchingGetBoostedResult>;

// Overload of different parameter type.
inline auto do_wbim_experiment_get_boosted(const WBIMAdjacencyListPair& graph, const VertexSet& seeds,
                                           const CommonExperimentParams& common, const WBIMSketchingParams& sketching)
    -> rfl::Result<exp_states::WBIMSketchingGetBoostedResult> {
  return exp_components::do_wbim_experiment_get_boosted( //
      graph, seeds, sketching, {.histogram_shape = common.histogram_shape()});
}

// ---- Group C: Simulation ----

// Preforms single simulation.
auto do_wim_simulate(const WIMAdjacencyListPair& graph, std::span<const vertex_id_t> seeds_selected,
                     uint64_t simulation_try_count) -> rfl::Result<exp_states::WIMSimulationInfo>;

// Performs multiple simulations.
// Let [K1, K2 ... Km] = n_seeds, the function gets simulation result with the first K1 seeds, the first K2 seeds, etc.
auto do_wim_simulate(const WIMAdjacencyListPair& graph, std::span<const vertex_id_t> seeds_selected,
                     std::span<const vertex_id_t> n_seeds, uint64_t simulation_try_count)
    -> rfl::Result<std::vector<exp_states::WIMSimulationInfo>>;

// Performs multiple simulations.
// Let [S1, S2 ... Sn] = seeds_selected, [K1, K2 ... Km] = n_seeds,
// The function gets simulation result with Si[0 .. Kj - 1] for each i = 1 to n, j = 1 to m.
// (n * m) simulations in total.
auto do_wim_experiment_simulate(const WIMAdjacencyListPair& graph, std::span<const VertexList> seeds_selected,
                                const WIMSketchingParams& params, uint64_t simulation_try_count)
    -> rfl::Result<std::vector<exp_states::WIMSketchingSimulationResult>>;

// Performs single simulation.
// Let B = the boosted set, base = the objective function when B is empty set.
// If base is known, let base_objective_function point to the value; Otherwise, let base_objective_function be nullptr.
auto do_wbim_simulate(const WBIMAdjacencyListPair& graph, const VertexSet& seeds,
                      std::span<const vertex_id_t> boosted_selected, uint64_t simulation_try_count,
                      const double* base_objective_function = nullptr) -> rfl::Result<exp_states::WBIMSimulationInfo>;

// Performs multiple simulations.
// Let [K1, K2 ... Km] = n_boosted,
// the function gets simulation result with the first K1 boosted vertices, the first K2 boosted vertices, etc.
// base_objective_function see above.
auto do_wbim_simulate(const WBIMAdjacencyListPair& graph, const VertexSet& seeds,
                      std::span<const vertex_id_t> boosted_selected, std::span<const vertex_id_t> n_boosted,
                      uint64_t simulation_try_count, const double* base_objective_function = nullptr)
    -> rfl::Result<std::vector<exp_states::WBIMSimulationInfo>>;

// Performs multiple simulations.
// Let [B1, B2 ... Bn] = boosted_selected, [K1, K2 ... Km] = n_boosted,
// The function gets simulation result with Bi[0 .. Kj - 1] for each i = 1 to n, j = 1 to m.
// (n * m) simulations in total.
auto do_wbim_experiment_simulate(const WBIMAdjacencyListPair& graph, const VertexSet& seeds,
                                 std::span<const VertexList> boosted_selected, const WBIMSketchingParams& params,
                                 uint64_t simulation_try_count, const double* base_objective_function = nullptr)
    -> rfl::Result<std::vector<exp_states::WBIMSketchingSimulationResult>>;

// ---- Group D: Graph Coarsening ----

// Continues coarsening, with newly coarsened levels appended to coarsen_results.
//
// If coarsen_results is empty, then coarsening starts from original_graph, whose level is 0.
// Otherwise, let L = coarsen_results.size(), coarsening starts from coarsen_results.back(), whose level is L.
//
// Coarsening stops until one of the following contidions is met:
//   1. Coarsening level reaches destination_level;
//   2. # of coarsened vertices is no more than coarsening_threshold.
//
// A dummy exp_states::CoarseningInfo object is returned (members all-zero)
// when no coarsening operation is performed.
// (i.e. original_graph or coarsened_results.back() already meets the conditions above)
auto do_continue_coarsening_wim_graph(const WIMAdjacencyListPair& original_graph,
                                      std::vector<WIMCoarsenGraphBriefResult>& coarsen_results,
                                      vertex_id_t destination_level, vertex_id_t coarsening_threshold,
                                      const CoarseningParams& params) -> rfl::Result<exp_states::CoarseningInfo>;

// Overload of different parameter type
inline auto do_continue_coarsening_wim_graph(const WIMAdjacencyListPair& original_graph,
                                             std::vector<WIMCoarsenGraphBriefResult>& coarsen_results,
                                             vertex_id_t destination_level,
                                             const MultiLevelParams& multi_level_params) {
  return exp_components::do_continue_coarsening_wim_graph( //
      original_graph, coarsen_results, destination_level,  //
      *multi_level_params.coarsening_threshold, *multi_level_params.coarsening);
}

// Same as above.
auto do_continue_coarsening_wbim_graph(const WBIMAdjacencyListPair& original_graph, const VertexSet& seeds,
                                       std::vector<WBIMCoarsenGraphBriefResult>& coarsen_results,
                                       vertex_id_t destination_level, vertex_id_t coarsening_threshold,
                                       const CoarseningParams& params) -> rfl::Result<exp_states::CoarseningInfo>;

// Overload of different parameter type
inline auto do_continue_coarsening_wbim_graph(const WBIMAdjacencyListPair& original_graph, const VertexSet& seeds,
                                              std::vector<WBIMCoarsenGraphBriefResult>& coarsen_results,
                                              vertex_id_t destination_level,
                                              const MultiLevelParams& multi_level_params) {
  return exp_components::do_continue_coarsening_wbim_graph(      //
      original_graph, seeds, coarsen_results, destination_level, //
      *multi_level_params.coarsening_threshold, *multi_level_params.coarsening);
}

// Coarsens the given graph until one of the following conditions is met:
//   1. Coarsening level reaches destination_level;
//   2. # of coarsened vertices is no more than coarsening_threshold.
auto do_coarsen_wim_graph(const WIMAdjacencyListPair& graph, vertex_id_t destination_level,
                          vertex_id_t coarsening_threshold, const CoarseningParams& params)
    -> rfl::Result<exp_states::CoarseningResult<WIMEdge>>;

// Overload of different parameter type.
inline auto do_coarsen_wim_graph(const WIMAdjacencyListPair& graph, vertex_id_t destination_level,
                                 const MultiLevelParams& multi_level_params) {
  return exp_components::do_coarsen_wim_graph( //
      graph, destination_level, *multi_level_params.coarsening_threshold, *multi_level_params.coarsening);
}

// Coarsens the given graph until one of the following conditions is met:
//   1. Coarsening level reaches destination_level;
//   2. # of coarsened vertices is no more than coarsening_threshold.
auto do_coarsen_wbim_graph(const WBIMAdjacencyListPair& graph, const VertexSet& seeds, vertex_id_t destination_level,
                           vertex_id_t coarsening_threshold, const CoarseningParams& params)
    -> rfl::Result<exp_states::CoarseningResult<WBIMEdge>>;

// Overload of different parameter type.
inline auto do_coarsen_wbim_graph(const WBIMAdjacencyListPair& graph, const VertexSet& seeds,
                                  vertex_id_t destination_level, const MultiLevelParams& multi_level_params) {
  return exp_components::do_coarsen_wbim_graph( //
      graph, seeds, destination_level, *multi_level_params.coarsening_threshold, *multi_level_params.coarsening);
}

// ---- Group E: Expanding

// Expands the given vertex list to original_graph.
// Let L = coarsen_results.size() be the total number of coarsening levels, l = n_fast_expanding_levels.
// The first l levels (0 -> 1 -> ... -> l) are expanded with O(K) fast method, where K = vertices.size().
// Other levels (l -> ... -> L) are expanded with given params.
auto do_expand_wim_vertices(const AdjacencyListPair<WIMEdge>& original_graph,
                            std::span<const WIMCoarsenGraphBriefResult> coarsen_results, VertexList vertices,
                            vertex_id_t n_fast_expanding_levels, const ExpandingParams& params)
    -> rfl::Result<exp_states::WIMExpansionInfo>;

// Expands each of the given vertex lists to original_graph.
// n_fast_expanding_levels and params see above.
auto do_expand_wim_vertex_lists(const AdjacencyListPair<WIMEdge>& original_graph,
                                std::span<const WIMCoarsenGraphBriefResult> coarsen_results,
                                VertexListList vertex_lists, vertex_id_t n_fast_expanding_levels,
                                const ExpandingParams& params) -> rfl::Result<exp_states::WIMExpansionResult>;

// Overload of different parameter type.
inline auto do_expand_wim_vertex_lists(const AdjacencyListPair<WIMEdge>& original_graph,
                                       std::span<const WIMCoarsenGraphBriefResult> coarsen_results,
                                       VertexListList vertex_lists, const MultiLevelParams& multi_level_params) {
  return exp_components::do_expand_wim_vertex_lists(            //
      original_graph, coarsen_results, std::move(vertex_lists), //
      *multi_level_params.n_fast_expanding_levels, *multi_level_params.expanding);
}

// Expands the given vertex list to original_graph.
// n_fast_expanding_levels and params see above.
auto do_expand_wbim_vertices(const AdjacencyListPair<WIMEdge>& original_graph,
                             std::span<const WIMCoarsenGraphBriefResult> coarsen_results, VertexList vertices,
                             vertex_id_t n_fast_expanding_levels, const ExpandingParams& params)
    -> rfl::Result<exp_states::WIMExpansionInfo>;

// Expands each of the given vertex lists to original_graph.
// n_fast_expanding_levels and params see above.
auto do_expand_wbim_vertex_lists(const AdjacencyListPair<WBIMEdge>& original_graph,
                                 std::span<const WBIMCoarsenGraphBriefResult> coarsen_results,
                                 VertexListList vertex_lists, vertex_id_t n_fast_expanding_levels,
                                 const ExpandingParams& params) -> rfl::Result<exp_states::WBIMExpansionResult>;

// Overload of different parameter type.
inline auto do_expand_wbim_vertex_lists(const AdjacencyListPair<WBIMEdge>& original_graph,
                                        std::span<const WBIMCoarsenGraphBriefResult> coarsen_results,
                                        VertexListList vertex_lists, const MultiLevelParams& multi_level_params) {
  return exp_components::do_expand_wbim_vertex_lists(           //
      original_graph, coarsen_results, std::move(vertex_lists), //
      *multi_level_params.n_fast_expanding_levels, *multi_level_params.expanding);
}
} // namespace exp_components
