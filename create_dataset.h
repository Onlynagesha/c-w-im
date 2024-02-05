#pragma once

#include "graph_types.h"
#include <rfl/Result.hpp>
#include <rfl/Validator.hpp>
#include <rfl/always_false.hpp>
#include <rfl/comparisons.hpp>
#include <ylt/easylog/record.hpp>

struct ReadGraphParams {
  std::string input_file;
  std::string output_file;

  rfl::Validator<double, rfl::Minimum<0>, rfl::Maximum<1>> mu_min = 0.0;
  rfl::Validator<double, rfl::Minimum<0>, rfl::Maximum<1>> mu_max = 1.0;

  rfl::Validator<double, rfl::Minimum<0>> lambda_in = 1.0;
  rfl::Validator<double, rfl::Minimum<0>> lambda_seed = 1.0;
  rfl::Validator<double, rfl::Minimum<0>> lambda_boost = 1.0;

  rfl::Validator<double, rfl::Minimum<1>> alpha = 1.0;
  rfl::Validator<double, rfl::Minimum<1>> alpha_seed = 1.0;
  rfl::Validator<double, rfl::Minimum<1>> alpha_boost = 1.0;

  // Vertices with larger degree will be forced to be unweighted.
  rfl::Validator<double, rfl::Minimum<0>, rfl::Maximum<1>> unweighted_ratio = 0.0;
  // Edges whose source vertex has larger degree will be penalized
  rfl::Validator<double, rfl::Minimum<0>, rfl::Maximum<1>> penalized_ratio = 0.0;
  // For each edge e = (u, v) whose source vertex u is penalized by the ratio above,
  // let p(e), p_seed(e), p_boost(e) /= penalize_factor.
  rfl::Validator<double, rfl::Minimum<1>> penalize_factor = 1.0;

  easylog::Severity log_level = easylog::Severity::INFO;
  rfl::Validator<size_t, rfl::Minimum<1>> histogram_width = 100;
  rfl::Validator<size_t, rfl::Minimum<1>> histogram_height = 20;

  static auto parse_from_args(int argc, char** argv) noexcept -> rfl::Result<ReadGraphParams>;

  auto get_wim_p_values(double in_deg, double out_deg, double penalize_factor = 1.0) const noexcept -> WIMEdge;
  auto get_wbim_p_values(double in_deg, double out_deg, double penalize_factor = 1.0) const noexcept -> WBIMEdge;

  template <is_edge_property E>
  auto get_p_values(double in_deg, double out_deg, double penalize_factor = 1.0) const noexcept {
    if constexpr (std::is_same_v<E, WIMEdge>) {
      return get_wim_p_values(in_deg, out_deg, penalize_factor);
    } else if constexpr (std::is_same_v<E, WBIMEdge>) {
      return get_wbim_p_values(in_deg, out_deg, penalize_factor);
    } else {
      static_assert(rfl::always_false_v<E>, "Invalid edge type.");
    }
  }
};

auto create_wim_dataset(const ReadGraphParams& params) noexcept -> WIMReadGraphResult;
auto create_wim_dataset_jstr(const std::string& params_json_text) noexcept -> rfl::Result<WIMReadGraphResult>;
