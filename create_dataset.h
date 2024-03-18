#pragma once

#include "graph_types.h"
#include <rfl/Flatten.hpp>
#include <rfl/Result.hpp>
#include <rfl/Validator.hpp>
#include <rfl/always_false.hpp>
#include <rfl/comparisons.hpp>
#include <ylt/easylog/record.hpp>

struct CreateDatasetCommonParams {
  std::string input_file;
  std::string output_file;

  rfl::Validator<double, rfl::Minimum<0>, rfl::Maximum<1>> mu_min = 0.0;
  rfl::Validator<double, rfl::Minimum<0>, rfl::Maximum<1>> mu_max = 1.0;

  rfl::Validator<double, rfl::Minimum<0>> lambda_in = 1.0;
  rfl::Validator<double, rfl::Minimum<1>> alpha = 1.0;

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
};

struct CreateWIMDatasetParams {
  rfl::Flatten<CreateDatasetCommonParams> common;
  rfl::Validator<double, rfl::Minimum<0>> lambda_seed = 1.0;
  rfl::Validator<double, rfl::Minimum<1>> alpha_seed = 1.0;

  static auto parse_from_args(int argc, char** argv) noexcept -> rfl::Result<CreateWIMDatasetParams>;

  auto make_edge(double in_deg, double out_deg, bool penalized) const noexcept -> WIMEdge;
};

struct CreateWBIMDatasetParams {
  rfl::Flatten<CreateDatasetCommonParams> common;
  rfl::Validator<double, rfl::Minimum<0>> lambda_boost = 1.0;
  rfl::Validator<double, rfl::Minimum<1>> alpha_boost = 1.0;

  static auto parse_from_args(int argc, char** argv) noexcept -> rfl::Result<CreateWBIMDatasetParams>;

  auto make_edge(double in_deg, double out_deg, bool penalized) const noexcept -> WBIMEdge;
};

template <is_edge_property E>
struct CreateDatasetParamsTraits {
  static_assert(rfl::always_false_v<E>, "Invalid edge type.");
};

template <>
struct CreateDatasetParamsTraits<WIMEdge> {
  using ParamsType = CreateWIMDatasetParams;
};

template <>
struct CreateDatasetParamsTraits<WBIMEdge> {
  using ParamsType = CreateWBIMDatasetParams;
};

auto create_wim_dataset(const CreateWIMDatasetParams& params) noexcept -> WIMReadGraphResult;
auto create_wim_dataset_jstr(const std::string& params_json_text) noexcept -> rfl::Result<WIMReadGraphResult>;

auto create_wbim_dataset(const CreateWBIMDatasetParams& params) noexcept -> WBIMReadGraphResult;
auto create_wbim_dataset_jstr(const std::string& params_json_text) noexcept -> rfl::Result<WBIMReadGraphResult>;

template <is_edge_property E>
inline auto create_dataset(const typename CreateDatasetParamsTraits<E>::ParamsType& params) noexcept {
  if constexpr (std::is_same_v<E, WIMEdge>) {
    return create_wim_dataset(params);
  } else if constexpr (std::is_same_v<E, WBIMEdge>) {
    return create_wbim_dataset(params);
  } else {
    static_assert(rfl::always_false_v<E>, "Invalid edge type.");
  }
}
