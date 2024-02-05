#include "dump.h"
#include "utils/result.h"
#include "wim.h"
#include <nlohmann/json.hpp>
#include <nwgraph/util/timer.hpp>
#include <ylt/easylog.hpp>

auto create_json_fout_ptr(const std::string& out_path) {
  auto json_fout_ptr = std::unique_ptr<std::ofstream>{};
  if (!out_path.empty()) {
    auto fout = std::ofstream{out_path};
    if (!fout.is_open()) {
      ELOGFMT(ERROR, "Failed to open JSON output file '{}'. Uses log output as fallback.", out_path);
    } else {
      ELOGFMT(INFO, "Successfully opens JSON output file '{}'.", out_path);
      json_fout_ptr = std::make_unique<std::ofstream>(std::move(fout));
    }
  } else {
    ELOG_INFO << "No JSON output file specified. Uses log output as fallback.";
  }
  return json_fout_ptr;
}

auto main_worker(int argc, char** argv) -> ResultVoid try {
  auto timer = nw::util::seconds_timer{};

  return WIMParams::parse_from_args(argc, argv).and_then([&](const WIMParams& params) {
    easylog::set_min_severity(params.log_severity);
    if (!params.log_output_file.empty()) {
      easylog::init_log(params.log_severity, params.log_output_file, false, params.log_console);
    }
    ELOGFMT(INFO, "Parameters: {:4}", params);
    auto json_fout_ptr = create_json_fout_ptr(params.json_output_file);

    timer.start();
    return read_wim_graph_data(params.input_file).and_then([&](const auto& read_result) {
      auto read_graph_time = timer.lap();
      ELOGFMT(INFO, "Done reading graph. |V| = {}, |E| = {}, time usage = {:.3} sec.",
              read_result.adj_list.num_vertices()[0], read_result.adj_list.num_edges(), read_graph_time);

      return wim_experiment(read_result, params).and_then([&](json json_root) -> ResultVoid {
        json_root["time_used"]["read_graph"] = read_graph_time;
        auto json_root_str = json_root.dump(4);
        if (json_fout_ptr) {
          ELOGFMT(DEBUG, "JSON output: {}", json_root_str);
          (*json_fout_ptr) << json_root_str;
        } else {
          ELOGFMT(INFO, "JSON output: {}", json_root_str);
        }
        return RESULT_VOID_SUCCESS;
      });
    });
  });
}
RFL_RESULT_CATCH_HANDLER()

int main(int argc, char** argv) {
  return main_worker(argc, argv)
      .and_then([](auto) {
        ELOG_INFO << "WIM experiment done.";
        return rfl::Result{0};
      })
      .or_else([](const rfl::Error& error) {
        ELOGFMT(CRITICAL, "WIM experiment error: `{}'", error.what());
        return rfl::Result{-1};
      })
      .value_or(-1);
}
