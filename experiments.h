#pragma once

#include "coarsening.h"
#include "contrast_algorithms.h"
#include "utils/histogram_shape.h"
#include "utils/result.h"
#include <rfl/Flatten.hpp>
#include <ylt/easylog.hpp>

// These parameters are expected to be "inherited" with rfl::Flatten
struct CommonExperimentParams {
  // Input file of graph data
  std::string input_file;
  // Output file of log.
  std::string log_output_file;
  // Output file of JSON which contains the experiment results.
  std::string json_output_file;
  // Log level, DEBUG or INFO recommended.
  easylog::Severity log_severity = easylog::Severity::DEBUG;
  // Whether or not to output log message to stdout.
  bool log_console = false;
  // Width of histograms during data distribution display.
  size_t histogram_width = 100;
  // Height of histograms during data distribution display.
  size_t histogram_height = 20;
  // # of Monte-carlo trials during simulation to estimate the objective function F(S) or F(B;S)
  rfl::Validator<uint64_t, rfl::Minimum<1>> simulation_try_count = 10'000;

  auto histogram_shape() const -> HistogramShape {
    return {.display_width = histogram_width, .display_height = histogram_height};
  }
};

// These parameters are expected to be "inherited" with rfl::Flatten
struct WIMSketchingParams {
  std::vector<size_t> n_sketches;
  std::vector<vertex_id_t> n_seeds;
};

struct WBIMSketchingParams {
  std::vector<size_t> n_sketches;
  std::vector<vertex_id_t> n_boosted;
};

struct WBIMSeedGeneratingParams {
  rfl::Rename<"n_seeds_to_generate", rfl::Validator<vertex_id_t, rfl::Minimum<1>>> //
      n_seeds = 10;
  rfl::Rename<"seed_candidate_ratio", rfl::Validator<double, rfl::ExclusiveMinimum<0>, rfl::Maximum<1>>> //
      candidate_ratio = 0.1;
  rfl::Validator<vertex_id_t, rfl::Minimum<1>> //
      max_n_wcc_without_seeds = 1;
  rfl::Rename<"seed_generation_max_try_count", rfl::Validator<uint64_t, rfl::Minimum<1>>> //
      max_try_count = 10'000;
};

struct WIMSketchingExperimentParams {
  rfl::Flatten<CommonExperimentParams> common;
  rfl::Flatten<WIMSketchingParams> sketching;

  static auto parse_from_args(int argc, char** argv) noexcept -> rfl::Result<WIMSketchingExperimentParams>;
};

struct WBIMSketchingExperimentParams {
  rfl::Flatten<CommonExperimentParams> common;
  rfl::Flatten<WBIMSketchingParams> sketching;
  rfl::Flatten<WBIMSeedGeneratingParams> seed_generating;

  static auto parse_from_args(int argc, char** argv) noexcept -> rfl::Result<WBIMSketchingExperimentParams>;
};

struct MultiLevelParams {
  rfl::Flatten<CoarseningParams> coarsening;
  rfl::Flatten<ExpandingParams> expanding;
  // Coarsening stops when |V| <= coarsening_threshold
  rfl::Validator<vertex_id_t, rfl::Minimum<1>> coarsening_threshold = 1'000;
  // LOCAL (which is the fastest) is the forced expanding policy at the first F levels
  // where F = n_fast_expanding_levels, since other policies are usually time-consuming.
  rfl::Validator<vertex_id_t, rfl::Minimum<0>> n_fast_expanding_levels = 0;
};

struct WIMCoarseningExperimentParams {
  rfl::Flatten<CommonExperimentParams> common;
  rfl::Flatten<WIMSketchingParams> sketching;
  rfl::Flatten<MultiLevelParams> multi_level;
  bool skips_first_level = false;

  static auto parse_from_args(int argc, char** argv) noexcept -> rfl::Result<WIMCoarseningExperimentParams>;
};

struct WBIMCoarseningExperimentParams {
  rfl::Flatten<CommonExperimentParams> common;
  rfl::Flatten<WBIMSketchingParams> sketching;
  rfl::Flatten<MultiLevelParams> multi_level;
  rfl::Flatten<WBIMSeedGeneratingParams> seed_generating;
  bool skips_first_level = false;

  static auto parse_from_args(int argc, char** argv) noexcept -> rfl::Result<WBIMCoarseningExperimentParams>;
};

struct WIMContrastExperimentParams {
  bool with_max_degree = false;
  bool with_max_strength = false;
  bool with_imrank = false;
  bool with_pagerank = false;
  bool with_rr_sketch = false;
  bool with_r_robust_scc = false;

  rfl::Flatten<CommonExperimentParams> common;
  std::vector<vertex_id_t> n_seeds;
  // For coarsening
  vertex_id_t coarsening_level = 0;
  rfl::Flatten<MultiLevelParams> multi_level;
  // For RR-sketching only
  std::vector<size_t> n_sketches;
  // For Pagerank only
  PagerankParams::DampingFactor pagerank_damping_factor = 0.85;
  PagerankParams::Epsilon pagerank_epsilon = 1e-6_ep;
  uint64_t pagerank_n_iterations = 100;
  // For IMRank only
  uint64_t imrank_n_iterations = 100;
  rfl::Validator<uint64_t, rfl::Minimum<1>> imrank_n_iterations_before_topk_fixed = 3;
  // For r-Robust SCC only
  uint64_t r_robust_scc_r = 16;

  static auto parse_from_args(int argc, char** argv) noexcept -> rfl::Result<WIMContrastExperimentParams>;
};

auto wim_experiment(int argc, char** argv) noexcept -> ResultVoid;

auto wbim_experiment(int argc, char** argv) noexcept -> ResultVoid;

auto wim_coarsening_experiment(int argc, char** argv) noexcept -> ResultVoid;

auto wbim_coarsening_experiment(int argc, char** argv) noexcept -> ResultVoid;

auto wim_contrast_experiment(int argc, char** argv) noexcept -> ResultVoid;
