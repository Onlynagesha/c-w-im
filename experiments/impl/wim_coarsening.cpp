#include "experiments/impl/wim_coarsening.h"
#include "../../wim.h"
#include "experiments/components.h"
#include "experiments/frameworks.h"
#include "experiments/impl/wim.h"

using exp_states::VertexList;
using exp_states::VertexListList;

namespace exp_impl {
auto do_wim_experiment_expand_simulate(const AdjacencyListPair<WIMEdge>& original_graph,
                                       std::span<const WIMCoarsenGraphBriefResult> coarsen_results,
                                       VertexListList coarsened_seed_lists, uint64_t simulation_try_count,
                                       const WIMSketchingParams& sketching_params,
                                       const MultiLevelParams& multi_level_params, json* json_root) -> ResultVoid {
  // Step 2: Expanding to original graph
  return exp_components::do_expand_wim_vertex_lists( //
             original_graph, coarsen_results, std::move(coarsened_seed_lists), multi_level_params)
      .and_then([&](exp_states::WIMExpansionResult expansion_res) {
        (*json_root)["expanding"] = exp_states::to_json(expansion_res);
        // Step 3: Simulation
        return exp_components::do_wim_experiment_simulate( //
            original_graph, expansion_res.expanded_seeds, sketching_params, simulation_try_count);
      })
      .transform([&](std::vector<exp_states::WIMSketchingSimulationResult> sim_res) {
        // Finally, outputs results to JSON.
        (*json_root)["simulating"] = exp_states::to_json(sim_res);
        return RESULT_VOID_SUCCESS;
      });
}

auto do_wim_experiment_with_coarsened(const AdjacencyListPair<WIMEdge>& original_graph,
                                      std::span<const WIMCoarsenGraphBriefResult> coarsen_results,
                                      const CommonExperimentParams& common_params,
                                      const WIMSketchingParams& sketching_params,
                                      const MultiLevelParams& multi_level_params, json* json_root) -> ResultVoid {
  BOOST_ASSERT(json_root != nullptr);
  const auto* graph_for_algo = coarsen_results.empty() ? &original_graph : &coarsen_results.back().coarsened;
  // Step 1: Get coarsened seeds by RR-sketching
  return exp_components::do_wim_experiment_get_seeds(*graph_for_algo, common_params, sketching_params)
      .and_then([&](exp_states::WIMSketchingGetSeedsResult coarsened_res) {
        (*json_root)["sketching"] = exp_states::to_json(coarsened_res);
        // Step 2 ~ 3: See above
        return do_wim_experiment_expand_simulate( //
            original_graph, coarsen_results, coarsened_res.selected_seeds, *common_params.simulation_try_count,
            sketching_params, multi_level_params, json_root);
      });
}

auto do_wbim_experiment_with_coarsened(const AdjacencyListPair<WBIMEdge>& original_graph, const VertexSet& seeds,
                                       std::span<const WBIMCoarsenGraphBriefResult> coarsen_results,
                                       const CommonExperimentParams& common_params,
                                       const WBIMSketchingParams& sketching_params,
                                       const MultiLevelParams& multi_level_params, json* json_root) -> ResultVoid {
  BOOST_ASSERT(json_root != nullptr);

  const auto* graph_for_sketching = //
      coarsen_results.empty() ? &original_graph : &coarsen_results.back().coarsened;
  const auto* seeds_for_sketching = //
      coarsen_results.empty() ? &seeds : &coarsen_results.back().coarsened_seeds;

  const auto& [orig_adj_list, orig_inv_adj_list, orig_vertex_weights] = original_graph;
  auto [orig_n, orig_m] = original_graph.graph_n_m();
  auto sim_try_count = *common_params.simulation_try_count;

  // Gets the base objective function when the boosted vertex set B is empty.
  auto base_F_obj = wbim_simulate_w(orig_adj_list, orig_vertex_weights, seeds, {orig_n, {}}, sim_try_count);
  RFL_RETURN_ON_ERROR(base_F_obj);
  auto base_F = *base_F_obj;

  // Step 1: Gets boosted vertices by PRR-sketching
  return exp_components::do_wbim_experiment_get_boosted( //
             *graph_for_sketching, *seeds_for_sketching, common_params, sketching_params)
      .and_then([&](exp_states::WBIMSketchingGetBoostedResult coarsened_res) -> ResultVoid {
        (*json_root)["sketching"] = exp_states::to_json(coarsened_res);

        auto do_step_2_3 = [&](std::string_view expanding_key, std::string_view simulating_key,
                               const VertexListList& selected_boosted) -> ResultVoid {
          // Step 2: Expands the boosted vertices to the original graph
          return exp_components::do_expand_wbim_vertex_lists( //
                     original_graph, coarsen_results, selected_boosted, multi_level_params)
              .and_then([&](exp_states::WBIMExpansionResult expansion_res) {
                (*json_root)[expanding_key] = exp_states::to_json(expansion_res);
                // Step 3: Simulation
                return exp_components::do_wbim_experiment_simulate( //
                    original_graph, seeds, expansion_res.expanded_boosted, sketching_params, sim_try_count, &base_F);
              })
              .transform([&](std::vector<exp_states::WBIMSketchingSimulationResult> sim_res) {
                // Finally, writes the result to JSON.
                (*json_root)[simulating_key] = exp_states::to_json(sim_res);
                return RESULT_VOID_SUCCESS;
              });
        };
        // Note: in WBIM experiments, two groups of boosted vertices are obtained:
        //   (1) via regular greedy selection algorithm;
        //   (2) via a faster algorithm which selects the best "critical vertices" as boosted.
        return do_step_2_3("expanding", "simulating", // (1)
                           coarsened_res.selected_boosted)
            .and_then([&](rfl::Nothing) {
              return do_step_2_3("expanding_by_critical", "simulating_by_critical", // (2)
                                 coarsened_res.selected_boosted_by_critical);
            });
      });
}

auto do_wim_coarsening_experiment(const AdjacencyListPair<WIMEdge>& graph, const WIMCoarseningExperimentParams& params,
                                  json* json_root) -> ResultVoid {
  BOOST_ASSERT(json_root != nullptr);
  auto timer = nw::util::seconds_timer{};
  auto& json_exp = (*json_root)["experiments"];

  const auto& [adj_list, inv_adj_list, vertex_weights] = graph;
  auto [n, m] = exp_io::write_graph_basic_information(*json_root, adj_list);

  if (n <= *params.multi_level->coarsening_threshold) {
    constexpr auto msg_pattern = "Experiments fails to proceed since |V| <= threshold, with |V| = {}, threshold = {}";
    return rfl::Error{fmt::format(msg_pattern, n, *params.multi_level->coarsening_threshold)};
  }
  // Step 1: Level 0, i.e. initial graph
  if (!params.skips_first_level) {
    ELOG_INFO << "Starts Level 0: solves WIM problem on the original graph.";
    auto& json_exp_cur = json_exp.emplace_back(json{{"level", 0}});
    auto connectivity_info = exp_io::write_graph_connectivity_information(json_exp_cur, graph);
    ELOGFMT(INFO, "# of weakly connected components = {}", connectivity_info.n_wcc);
    RFL_RETURN_ON_ERROR(do_wim_experiment(graph, *params.common, *params.sketching, &json_exp_cur));
  }

  // Step 2: Coarsening until the size threshold
  auto coarsen_results = std::vector<WIMCoarsenGraphBriefResult>{};
  for (auto level = 1_vid, n_coarsened = n; n_coarsened > params.multi_level->coarsening_threshold; ++level) {
    auto& json_exp_cur = json_exp.emplace_back(json{{"level", level}});
    // Step 2.1: Gets the coarsened graph
    auto coarsen_success = //
        exp_components::do_continue_coarsening_wim_graph(graph, coarsen_results, level, *params.multi_level)
            .transform([&](exp_states::CoarseningInfo info) {
              json_exp_cur["coarsening"] = exp_states::to_json(info);
              return RESULT_VOID_SUCCESS;
            });
    RFL_RETURN_ON_ERROR(coarsen_success);
    const auto& cur_coarsen_result = coarsen_results.back();
    n_coarsened = cur_coarsen_result.coarsened.n_vertices();

    // Checks connectivity for debugging & data analysis
    auto connectivity_info = exp_io::write_graph_connectivity_information(json_exp_cur, cur_coarsen_result.coarsened);
    ELOGFMT(INFO, "# of weakly connected components = {}", connectivity_info.n_wcc);
    // Step 2.2: Solves on the coarsened graph -> exp_states::anding seeds -> Simulation to estimate F(S)
    RFL_RETURN_ON_ERROR(do_wim_experiment_with_coarsened(graph, coarsen_results, params, &json_exp_cur));
  }
  return RESULT_VOID_SUCCESS;
}

auto do_wbim_coarsening_experiment(const AdjacencyListPair<WBIMEdge>& graph, const VertexSet& seeds,
                                   const WBIMCoarseningExperimentParams& params, json* json_root) -> ResultVoid {
  BOOST_ASSERT(json_root != nullptr);
  auto timer = nw::util::seconds_timer{};
  auto& json_exp = (*json_root)["experiments"];

  const auto& [adj_list, inv_adj_list, vertex_weights] = graph;
  auto [n, m] = exp_io::write_graph_basic_information(*json_root, adj_list);

  if (n <= *params.multi_level->coarsening_threshold) {
    constexpr auto msg_pattern = "Experiments fails to proceed since |V| <= threshold, with |V| = {}, threshold = {}";
    return rfl::Error{fmt::format(msg_pattern, n, *params.multi_level->coarsening_threshold)};
  }
  // Step 1: Level 0, i.e. initial graph
  if (!params.skips_first_level) {
    ELOG_INFO << "Starts Level 0: solves WBIM problem on the original graph.";
    auto& json_exp_cur = json_exp.emplace_back(json{{"level", 0}});
    auto connectivity_info = exp_io::write_graph_connectivity_information(json_exp_cur, graph);
    ELOGFMT(INFO, "# of weakly connected components = {}", connectivity_info.n_wcc);
    RFL_RETURN_ON_ERROR(do_wbim_experiment(graph, seeds, *params.common, *params.sketching, &json_exp_cur));
  }

  // Step 2: Coarsening until the size threshold
  auto coarsen_results = std::vector<WBIMCoarsenGraphBriefResult>{};
  for (auto level = 1_vid, n_coarsened = n; n_coarsened > params.multi_level->coarsening_threshold; ++level) {
    auto& json_exp_cur = json_exp.emplace_back(json{{"level", level}});
    // Step 2.1: Gets the coarsened graph
    auto coarsen_success = //
        exp_components::do_continue_coarsening_wbim_graph(graph, seeds, coarsen_results, level, *params.multi_level)
            .transform([&](exp_states::CoarseningInfo info) {
              json_exp_cur["coarsening"] = exp_states::to_json(info);
              return RESULT_VOID_SUCCESS;
            });
    RFL_RETURN_ON_ERROR(coarsen_success);
    const auto& cur_coarsen_result = coarsen_results.back();
    n_coarsened = cur_coarsen_result.coarsened.n_vertices();

    // Checks connectivity for debugging & data analysis
    auto connectivity_info = exp_io::write_graph_connectivity_information(json_exp_cur, cur_coarsen_result.coarsened);
    ELOGFMT(INFO, "# of weakly connected components = {}", connectivity_info.n_wcc);
    // Step 2.2: Solves on the coarsened graph -> exp_states::anding seeds -> Simulation to estimate F(S)
    RFL_RETURN_ON_ERROR(do_wbim_experiment_with_coarsened(graph, seeds, coarsen_results, params, &json_exp_cur));
  }
  return RESULT_VOID_SUCCESS;
}
} // namespace exp_impl

auto wim_coarsening_experiment(int argc, char** argv) noexcept -> ResultVoid {
  return exp_frameworks::wim_experiment_framework<WIMCoarseningExperimentParams>( //
      argc, argv, exp_impl::do_wim_coarsening_experiment);
}

auto wbim_coarsening_experiment(int argc, char** argv) noexcept -> ResultVoid {
  return exp_frameworks::wbim_experiment_framework<WBIMCoarseningExperimentParams>( //
      argc, argv, exp_impl::do_wbim_coarsening_experiment);
}
