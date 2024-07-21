#pragma once

#include "dump.h"
#include "experiments.h"
#include "experiments/components.h"
#include "experiments/io.h"

namespace exp_frameworks {
template <class ParamsType, class DoExperimentFn>
  requires(std::is_invocable_r_v<ResultVoid, DoExperimentFn, // do_experiment_fn(graph, params, &json_root)
                                 const WIMAdjacencyListPair&, const ParamsType&, json*>)
auto wim_experiment_framework(int argc, char** argv, DoExperimentFn&& do_experiment_fn) -> ResultVoid try {
  auto timer = nw::util::seconds_timer{};

  // Step 1: Parses arguments from (argc, argv)
  return ParamsType::parse_from_args(argc, argv).and_then([&](ParamsType params) {
    exp_io::init_easylog(*params.common);
    ELOGFMT(INFO, "Parameters: {:4}", params);
    auto json_fout_ptr = exp_io::create_fout_ptr(params.common->json_output_file);

    timer.start();
    // Step 2: Reads graph data from input file
    return read_directed_wim_adjacency_lists(params.common->input_file)
        .and_then([&](AdjacencyListPair<WIMEdge> read_result) {
          auto read_graph_time = timer.lap();
          auto [n, m] = read_result.graph_n_m();
          ELOGFMT(INFO, "Done reading graph. |V| = {}, |E| = {}, time usage = {:.3} sec.", n, m, read_graph_time);

          auto json_root = json{};
          // Step 3: Performs experiment by calling do_experiment_fn
          return std::invoke(do_experiment_fn, read_result, params, &json_root).transform([&](rfl::Nothing) {
            auto json_root_str = json_root.dump(4);
            // Step 4: Outputs results as JSON to specified file
            exp_io::dump_to_fout_ptr(json_fout_ptr.get(), json_root_str);
            return RESULT_VOID_SUCCESS;
          });
        });
  });
}
RFL_RESULT_CATCH_HANDLER() // Error handling on exceptions

template <class ParamsType, class DoExperimentFn>
  requires(std::is_invocable_r_v<ResultVoid, DoExperimentFn, // do_experiment_fn(graph, seeds, params, &json_root)
                                 const WBIMAdjacencyListPair&, const VertexSet&, const ParamsType&, json*>)
auto wbim_experiment_framework(int argc, char** argv, DoExperimentFn&& do_experiment_fn) -> ResultVoid try {
  auto timer = nw::util::seconds_timer{};

  // Step 1: Parses arguments from (argc, argv)
  return ParamsType::parse_from_args(argc, argv).and_then([&](ParamsType params) {
    exp_io::init_easylog(*params.common);
    ELOGFMT(INFO, "Parameters: {:4}", params);
    auto json_fout_ptr = exp_io::create_fout_ptr(params.common->json_output_file);

    timer.start();
    // Step 2: Reads graph data from input file
    return read_directed_wbim_adjacency_lists(params.common->input_file)
        .and_then([&](AdjacencyListPair<WBIMEdge> read_result) -> ResultVoid {
          auto [n, m] = read_result.graph_n_m();
          // Step 3: Generates seed vertices for WBIM experiments
          return exp_components::generate_wbim_seeds(read_result, *params.seed_generating)
              .and_then([&](VertexSet seeds) {
                timer.stop();
                auto read_graph_time = timer.elapsed();
                constexpr auto msg_pattern = "Done reading graph and generating seeds for WBIM Experiment. "
                                             "|V| = {}, |E| = {}, time usage = {:.3} sec.";
                ELOGFMT(INFO, msg_pattern, n, m, read_graph_time);

                auto json_root = json{{"seeds", seeds.vertex_list}};
                // Step 4: Performs experiment by calling do_experiment_fn
                return std::invoke(do_experiment_fn, read_result, seeds, params, &json_root)
                    .transform([&](rfl::Nothing) {
                      auto json_root_str = json_root.dump(4);
                      // Step 5: Outputs results as JSON to specified file
                      exp_io::dump_to_fout_ptr(json_fout_ptr.get(), json_root_str);
                      return RESULT_VOID_SUCCESS;
                    });
              });
        });
  });
}
RFL_RESULT_CATCH_HANDLER()
} // namespace exp_frameworks