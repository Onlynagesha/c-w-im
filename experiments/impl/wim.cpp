#include "experiments/impl/wim.h"
#include "../../wim.h"
#include "experiments.h"
#include "experiments/components.h"
#include "experiments/frameworks.h"
#include "experiments/io.h"

namespace exp_impl {
auto do_wim_experiment(const WIMAdjacencyListPair& graph, const CommonExperimentParams& common,
                       const WIMSketchingParams& sketching, json* json_root) -> ResultVoid {
  BOOST_ASSERT(json_root != nullptr);
  exp_io::write_graph_basic_information(*json_root, graph.adj_list);

  // Step 1: Get seeds by RR-sketching
  return exp_components::do_wim_experiment_get_seeds(graph, common, sketching)
      .and_then([&](exp_states::WIMSketchingGetSeedsResult get_seeds_res) {
        // Step 2: Simulation
        return exp_components::do_wim_experiment_simulate( //
                   graph, get_seeds_res.selected_seeds, sketching, *common.simulation_try_count)
            .transform([&](std::vector<exp_states::WIMSketchingSimulationResult> sim_res) {
              BOOST_ASSERT(sim_res.size() == sketching.n_sketches.size());
              // Finally, outputs results to JSON.
              (*json_root)["sketching"] = exp_states::to_json(get_seeds_res);
              (*json_root)["simulating"] = exp_states::to_json(sim_res);
              return RESULT_VOID_SUCCESS;
            });
      });
}

auto do_wbim_experiment(const WBIMAdjacencyListPair& graph, const VertexSet& seeds,
                        const CommonExperimentParams& common, const WBIMSketchingParams& sketching, json* json_root)
    -> ResultVoid {
  BOOST_ASSERT(json_root != nullptr);
  exp_io::write_graph_basic_information(*json_root, graph.adj_list);

  // Step 1: Gets boosted vertices by PRR-sketching
  return exp_components::do_wbim_experiment_get_boosted(graph, seeds, common, sketching)
      .and_then([&](exp_states::WBIMSketchingGetBoostedResult get_boosted_res) {
        (*json_root)["sketching"] = exp_states::to_json(get_boosted_res);

        auto [n, m] = graph.graph_n_m();
        auto try_count = *common.simulation_try_count;
        // Step 2: Simulation
        return wbim_simulate_w(graph.adj_list, graph.vertex_weights, seeds, {n, {}}, try_count)
            .and_then([&](double base_F) -> ResultVoid {
              // Note: in WBIM experiments, two groups of boosted vertices are obtained:
              //   (1) via regular greedy selection algorithm;
              //   (2) via a faster algorithm which selects the best "critical vertices" as boosted.
              auto sim_res = exp_components::do_wbim_experiment_simulate( //
                  graph, seeds, get_boosted_res.selected_boosted, sketching, try_count, &base_F);
              RFL_RETURN_ON_ERROR(sim_res);
              auto sim_res_by_critical = exp_components::do_wbim_experiment_simulate( //
                  graph, seeds, get_boosted_res.selected_boosted_by_critical, sketching, try_count, &base_F);
              RFL_RETURN_ON_ERROR(sim_res_by_critical);
              // Finally, outputs results to JSON.
              (*json_root)["simulating"] = exp_states::to_json(*sim_res);
              (*json_root)["simulating_by_critical"] = exp_states::to_json(*sim_res_by_critical);
              return RESULT_VOID_SUCCESS;
            });
      });
}
} // namespace exp_impl

auto wim_experiment(int argc, char** argv) noexcept -> ResultVoid {
  auto do_experiment_fn = [](const AdjacencyListPair<WIMEdge>& graph, //
                             const WIMSketchingExperimentParams& params, json* json_root) {
    return exp_impl::do_wim_experiment(graph, *params.common, *params.sketching, json_root);
  };
  return exp_frameworks::wim_experiment_framework<WIMSketchingExperimentParams>(argc, argv, do_experiment_fn);
}

auto wbim_experiment(int argc, char** argv) noexcept -> ResultVoid {
  auto do_experiment_fn = [](const AdjacencyListPair<WBIMEdge>& graph, const VertexSet& seeds,
                             const WBIMSketchingExperimentParams& params, json* json_root) {
    return exp_impl::do_wbim_experiment(graph, seeds, *params.common, *params.sketching, json_root);
  };
  return exp_frameworks::wbim_experiment_framework<WBIMSketchingExperimentParams>(argc, argv, do_experiment_fn);
}
