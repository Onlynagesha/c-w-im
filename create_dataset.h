#pragma once

#include "graph_types.h"
#include "utils/boost_assert.h"
#include "utils/utils.h"
#include <rfl/Result.hpp>
#include <rfl/Validator.hpp>
#include <rfl/always_false.hpp>
#include <rfl/comparisons.hpp>

struct ReadGraphParams {
  std::string input_file;
  std::string output_file;

  rfl::Validator<double, rfl::Minimum<0>> lambda_in = 1.0;
  rfl::Validator<double, rfl::Minimum<0>> lambda_out = 1.0;

  rfl::Validator<double, rfl::Minimum<1>> alpha = 1.0;
  rfl::Validator<double, rfl::Minimum<1>> alpha_seed = 1.0;
  rfl::Validator<double, rfl::Minimum<1>> alpha_boost = 1.0;

  rfl::Validator<double, rfl::Minimum<0>, rfl::Maximum<1>> weighted_ratio;

  static auto parse_from_args(int argc, char** argv) noexcept -> rfl::Result<ReadGraphParams>;

  template <is_edge_property E>
  auto get_p_values(edge_probability_t in_deg, edge_probability_t out_deg) const {
    BOOST_ASSERT_MSG(in_deg > 0, "In-degree must be a positive integer.");
    BOOST_ASSERT_MSG(out_deg > 0, "Out-degree must be a positive integer.");

    auto p0 = std::min(lambda_in.value() / in_deg + lambda_out.value() / out_deg, 1.0);
    auto [p, p_seed, p_boost] = tuple_transform(std::tie(alpha, alpha_seed, alpha_boost), [&](auto a) {
      return (edge_probability_t)(1.0 - std::pow(1.0 - p0, a.value()));
    });
    if constexpr (std::same_as<E, WIMEdge>) {
      return WIMEdge{.p = p, .p_seed = p_seed};
    } else if constexpr (std::same_as<E, WBIMEdge>) {
      return WBIMEdge{.p = p, .p_seed = p_seed, .p_boost = p_boost};
    } else {
      static_assert(rfl::always_false_v<E>, "Invalid edge property type.");
    }
  }
};

auto create_wim_dataset(const ReadGraphParams& params) noexcept -> WIMReadGraphResult;
auto create_wim_dataset_jstr(const std::string& params_json_text) noexcept -> rfl::Result<WIMReadGraphResult>;
