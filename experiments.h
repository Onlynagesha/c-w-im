#pragma once

#include "coarsening.h"
#include "utils/histogram_shape.h"
#include "utils/result.h"
#include <rfl/Flatten.hpp>
#include <ylt/easylog.hpp>

// These parameters are expected to be "inherited" with rfl::Flatten
struct CommonExperimentParams {
  std::string input_file; // Input file of graph data
  std::string log_output_file;
  std::string json_output_file;
  easylog::Severity log_severity = easylog::Severity::DEBUG;
  bool log_console = false;
  size_t histogram_width = 100;
  size_t histogram_height = 20;

  auto histogram_shape() const -> HistogramShape {
    return {.display_width = histogram_width, .display_height = histogram_height};
  }
};

// These parameters are expected to be "inherited" with rfl::Flatten
struct WIMParams {
  std::vector<size_t> num_sketches;
  std::vector<vertex_id_t> num_seeds;
  rfl::Validator<uint64_t, rfl::Minimum<1>> simulation_try_count = 10'000;
};

struct WIMExperimentParams {
  rfl::Flatten<CommonExperimentParams> common;
  rfl::Flatten<WIMParams> wim;

  static auto parse_from_args(int argc, char** argv) noexcept -> rfl::Result<WIMExperimentParams>;
};

struct WIMCoarseningExperimentParams {
  rfl::Flatten<CommonExperimentParams> common;
  rfl::Flatten<WIMParams> wim;
  rfl::Flatten<CoarseningParams> coarsening;
  rfl::Flatten<ExpandingParams> expanding;
  // Coarsening stops when |V| <= coarsening_threshold
  rfl::Validator<vertex_id_t, rfl::Minimum<1>> coarsening_threshold = 1'000;
  // S_LOCAL (which is the fastest) is forced as expanding policy at the first F levels
  // with F = n_fast_expanding_levels
  rfl::Validator<vertex_id_t, rfl::Minimum<0>> n_fast_expanding_levels = 0;

  static auto parse_from_args(int argc, char** argv) noexcept -> rfl::Result<WIMCoarseningExperimentParams>;
};

auto wim_experiment(int argc, char** argv) noexcept -> ResultVoid;

auto wim_coarsening_experiment(int argc, char** argv) noexcept -> ResultVoid;
