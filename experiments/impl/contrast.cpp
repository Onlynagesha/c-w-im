#include "experiments.h"
#include "experiments/frameworks.h"
#include "experiments/impl/wim_coarsening.h"
#include <fmt/ranges.h>

namespace {
using exp_states::VertexList;
using exp_states::VertexListList;

auto do_wim_contrast_experiment(const AdjacencyListPair<WIMEdge>& graph, const WIMContrastExperimentParams& params,
                                json* json_root) -> ResultVoid {
  auto timer = nw::util::seconds_timer{};
  auto [n, m] = exp_io::write_graph_basic_information(*json_root, graph.adj_list);

  auto coarsen_results = //
      exp_components::do_coarsen_wim_graph(graph, params.coarsening_level, *params.multi_level)
          .transform([&](exp_states::CoarseningResult<WIMEdge> results) {
            (*json_root)["coarsening"] = exp_states::to_json(results.info);
            return std::move(results.coarsen_results);
          });
  RFL_RETURN_ON_ERROR(coarsen_results);

  if (params.with_rr_sketch) {
    auto& json_exp_cur = (*json_root)["rr_sketching"];
    auto sketching_params = WIMSketchingParams{
        .n_sketches = params.n_sketches,
        .n_seeds = params.n_seeds,
    };
    ELOG_INFO << "==== Starts RR-sketching algorithm ====";
    RFL_RETURN_ON_ERROR(exp_impl::do_wim_experiment_with_coarsened( //
        graph, *coarsen_results, *params.common, sketching_params, *params.multi_level, &json_exp_cur));
  }

#define CONTRAST_EXPERIMENT_FN [&](const WIMAdjacencyListPair& graph) -> rfl::Result<VertexList>

  auto contrast_experiment_framework = [&](std::string_view algorithm_name, auto&& algorithm_fn) -> ResultVoid {
    ELOGFMT(INFO, "==== Starts contrast algorithm: {} ====", algorithm_name);
    auto& json_exp_cur = (*json_root)[algorithm_name];
    timer.start();
    // Step 1: Performs the algorithm
    const auto* graph_for_algo = coarsen_results->empty() ? &graph : &coarsen_results->back().coarsened;
    return std::invoke(algorithm_fn, *graph_for_algo).and_then([&](VertexList seeds) {
      timer.stop();
      auto algo_time_used = timer.elapsed();
      json_exp_cur["algorithm"] = {
          {"seeds_before_expanding", seeds},
          {"time_used", algo_time_used},
      };
      ELOGFMT(INFO, "Finished algorithm '{}' in {:.3f} sec.", algorithm_name, algo_time_used);
      MYLOG_FMT_DEBUG("Seed vertices before expansion = {}", seeds);
      // Step 2 ~ 3: Expanding & Simulating
      auto sketching_params = WIMSketchingParams{
          .n_sketches = params.n_sketches,
          .n_seeds = params.n_seeds,
      };
      return exp_impl::do_wim_experiment_expand_simulate( //
          graph, *coarsen_results, {seeds},               //
          *params.common->simulation_try_count, sketching_params, *params.multi_level, &json_exp_cur);
    });
  };

  auto max_n_seeds = ranges::max(params.n_seeds);
  if (params.with_max_degree) {
    RFL_RETURN_ON_ERROR(contrast_experiment_framework(
        "max_degree", CONTRAST_EXPERIMENT_FN { return max_out_degree(graph.adj_list, max_n_seeds); }));
  }
  if (params.with_max_strength) {
    RFL_RETURN_ON_ERROR(contrast_experiment_framework(
        "max_strength", CONTRAST_EXPERIMENT_FN { return max_out_strength(graph.adj_list, max_n_seeds); }));
  }
  if (params.with_pagerank) {
    auto pagerank_params = PagerankParams{
        .damping_factor = *params.pagerank_damping_factor,
        .epsilon = *params.pagerank_epsilon,
        .k = max_n_seeds,
        .n_iterations = params.pagerank_n_iterations,
        .transpose = true,
        .uses_vertex_weight = true,
        .uses_edge_weight = true,
    };
    RFL_RETURN_ON_ERROR(contrast_experiment_framework(
        "pagerank", CONTRAST_EXPERIMENT_FN { return max_pagerank(graph, pagerank_params); }));
  }
  if (params.with_imrank) {
    auto imrank_params = IMRankParams{
        .k = max_n_seeds,
        .n_iterations = params.imrank_n_iterations,
        .n_iterations_before_topk_fixed = *params.imrank_n_iterations_before_topk_fixed,
    };
    RFL_RETURN_ON_ERROR(contrast_experiment_framework(
        "imrank", CONTRAST_EXPERIMENT_FN { return imrank(graph.inv_adj_list, imrank_params); }));
  }
  if (params.with_r_robust_scc) {
    auto r_robust_scc_params = RRobustSCCParams{
        .k = max_n_seeds,
        .r = params.r_robust_scc_r,
        .n_sketches = ranges::max(params.n_sketches),
    };
    RFL_RETURN_ON_ERROR(contrast_experiment_framework(
        "r_robust_scc", CONTRAST_EXPERIMENT_FN { return wim_r_robust_scc_p(graph, r_robust_scc_params); }));
  }
  // Done all
  return RESULT_VOID_SUCCESS;
}

#undef CONTRAST_EXPERIMENT_FN
} // namespace

auto wim_contrast_experiment(int argc, char** argv) noexcept -> ResultVoid {
  return exp_frameworks::wim_experiment_framework<WIMContrastExperimentParams>( //
      argc, argv, do_wim_contrast_experiment);
}
