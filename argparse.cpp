#include "create_dataset.h"
#include "experiments.h"
#include "to_matrix_market.h"
#include "utils/argparse_helper.h"
#include "utils/reflect_for_each_field.h"
#include "wim.h"
#include <fstream>

namespace {
struct ManualConfigArguments {
  bool has_config = false;
  bool requires_input_file = false;
  bool requires_output_file = false;
};

template <ManualConfigArguments Config>
auto manual_add_config_arguments(ArgumentParser& parser) -> void {
  if constexpr (Config.has_config) {
    parser.add_argument("-c", "--config").help("Configuration file as JSON format");
  }
  if constexpr (Config.requires_input_file) {
    parser.at("--input-file").required(); // Exception may be thrown
  }
  if constexpr (Config.requires_output_file) {
    parser.at("--output-file").required(); // Exception may be thrown
  }
}

template <class T>
auto read_arguments_from_config_file(ArgumentParser& parser, T& value) {
  if (auto in_path = parser.present("--config")) {
    read_fields_from_json(value, json::parse(std::ifstream{*in_path}));
  }
}

auto manual_sort_and_check_wim_params(WIMParams& params) -> std::optional<std::string> {
  auto do_sort_and_check = [](std::string_view name, auto& values) -> std::optional<std::string> {
    if (values.empty()) {
      return fmt::format("{} can not be an empty list.", name);
    }
    ranges::sort(values);
    if (auto n0 = values.front(); n0 <= 0) {
      return fmt::format("{} must be positive integers, while {} is given.", name, n0);
    }
    constexpr auto adjacent_equal_filter = views::filter(LAMBDA_1(get<0>(_1) == get<1>(_1)));
    for (auto [n0, n1] : views::adjacent<2>(values) | adjacent_equal_filter) {
      return fmt::format("Duplicated {} (detected {}) is disallowed.", name, n0);
    }
    return std::nullopt;
  };
  return do_sort_and_check("# of RR-sketches", params.num_sketches)
      .or_else(LAMBDA_0(do_sort_and_check("# of seed vertices", params.num_seeds)));
}
} // namespace

auto ReadGraphParams::parse_from_args(int argc, char** argv) noexcept -> rfl::Result<ReadGraphParams> {
  return parse_from_args_generic<ReadGraphParams>(
             argc, argv, // Appends --config etc.
             manual_add_config_arguments<ManualConfigArguments{
                 .has_config = true, .requires_input_file = true, .requires_output_file = true}>,
             read_arguments_from_config_file<ReadGraphParams>)
      .and_then([](ReadGraphParams params) -> rfl::Result<ReadGraphParams> {
        if (*params.lambda_in <= 0 && *params.lambda_seed <= 0) {
          return rfl::Error{"At least one of lambda_in and lambdas_seed shall be positive."};
        }
        if (*params.lambda_in <= 0 && *params.lambda_boost <= 0) {
          return rfl::Error{"At least one of lambda_in and lambda_boost shall be positive."};
        }
        return std::move(params);
      });
}

auto ToMatrixMarketParams::parse_from_args(int argc, char** argv) noexcept -> rfl::Result<ToMatrixMarketParams> {
  return parse_from_args_generic<ToMatrixMarketParams>(
             argc, argv, // Appends --config etc.
             manual_add_config_arguments<ManualConfigArguments{
                 .has_config = true, .requires_input_file = true, .requires_output_file = true}>,
             read_arguments_from_config_file<ToMatrixMarketParams>)
      .and_then([](ToMatrixMarketParams params) -> rfl::Result<ToMatrixMarketParams> {
        if (params.is_directed && params.is_undirected) {
          return rfl::Error{"--is-directed and --is-undirected shall not be set at the same time."};
        }
        return std::move(params);
      });
}

auto WIMExperimentParams::parse_from_args(int argc, char** argv) noexcept -> rfl::Result<WIMExperimentParams> {
  return parse_from_args_generic<WIMExperimentParams>(
             argc, argv, // Appends --config etc.
             manual_add_config_arguments<ManualConfigArguments{.has_config = true, .requires_input_file = true}>,
             read_arguments_from_config_file<WIMExperimentParams>)
      .and_then([](WIMExperimentParams params) -> rfl::Result<WIMExperimentParams> {
        if (auto err_msg = manual_sort_and_check_wim_params(*params.wim); err_msg) {
          return rfl::Error{std::move(*err_msg)};
        }
        return std::move(params);
      });
}

auto WIMCoarseningExperimentParams::parse_from_args(int argc, char** argv) noexcept
    -> rfl::Result<WIMCoarseningExperimentParams> {
  return parse_from_args_generic<WIMCoarseningExperimentParams>(
             argc, argv, // Appends --config etc.
             manual_add_config_arguments<ManualConfigArguments{.has_config = true, .requires_input_file = true}>,
             read_arguments_from_config_file<WIMCoarseningExperimentParams>)
      .and_then([](WIMCoarseningExperimentParams params) -> rfl::Result<WIMCoarseningExperimentParams> {
        if (auto err_msg = manual_sort_and_check_wim_params(*params.wim); err_msg) {
          return rfl::Error{std::move(*err_msg)};
        }
        return std::move(params);
      });
}
