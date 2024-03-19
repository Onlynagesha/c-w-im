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

struct CheckListArguments {
  bool non_empty_required = false;
  bool uniqueness_required = false;
  bool positive_values_required = false;
};

template <CheckListArguments Config, class T>
auto sort_and_check_list(std::string_view name, std::span<T> values) -> std::optional<std::string> {
  if (values.empty()) {
    if constexpr (Config.non_empty_required) {
      return fmt::format("'{}' can not be an empty list.", name);
    } else {
      return std::nullopt; // OK, no error message
    }
  }
  ranges::sort(values);
  if constexpr (Config.positive_values_required) {
    if (auto n0 = values.front(); n0 <= 0) {
      return fmt::format("Elements in '{}' must be positive, while {} is given.", name, n0);
    }
  }
  if constexpr (Config.uniqueness_required) {
    constexpr auto adjacent_equal_filter = views::filter(LAMBDA_1(get<0>(_1) == get<1>(_1)));
    for (auto [n0, n1] : views::adjacent<2>(values) | adjacent_equal_filter) {
      return fmt::format("Duplicated elements in '{}' (detected {}) is disallowed.", name, n0);
    }
  }
  return std::nullopt; // OK, no error
}

template <class SketchingParams>
auto sort_and_check_sketching_params(SketchingParams& sketching) -> std::optional<std::string> {
  constexpr auto requirements = CheckListArguments{
      .non_empty_required = true,
      .uniqueness_required = true,
      .positive_values_required = true,
  };
  auto n_sketches = std::span{sketching.n_sketches};
  return sort_and_check_list<requirements>("# of RR-sketches", n_sketches).or_else([&]() {
    if constexpr (std::is_same_v<SketchingParams, WIMSketchingParams>) {
      // n_seeds for WIM
      auto n_seeds = std::span{sketching.n_seeds};
      return sort_and_check_list<requirements>("# of seed vertices", n_seeds);
    } else if constexpr (std::is_same_v<SketchingParams, WBIMSketchingParams>) {
      // n_boosted for WBIM
      auto n_boosted = std::span{sketching.n_boosted};
      return sort_and_check_list<requirements>("# of boosted vertices", n_boosted);
    } else {
      static_assert(rfl::always_false_v<SketchingParams>, "Invalid sketching params type.");
    }
  });
}
} // namespace

auto CreateWIMDatasetParams::parse_from_args(int argc, char** argv) noexcept -> rfl::Result<CreateWIMDatasetParams> {
  return parse_from_args_generic<CreateWIMDatasetParams>(
             argc, argv, // Appends --config etc.
             manual_add_config_arguments<ManualConfigArguments{
                 .has_config = true, .requires_input_file = true, .requires_output_file = true}>,
             read_arguments_from_config_file<CreateWIMDatasetParams>)
      .and_then([](CreateWIMDatasetParams params) -> rfl::Result<CreateWIMDatasetParams> {
        if (*params.common->lambda_in <= 0 && *params.lambda_seed <= 0) {
          return rfl::Error{"At least one of lambda_in and lambdas_seed shall be positive."};
        }
        return std::move(params);
      });
}

auto CreateWBIMDatasetParams::parse_from_args(int argc, char** argv) noexcept -> rfl::Result<CreateWBIMDatasetParams> {
  return parse_from_args_generic<CreateWBIMDatasetParams>(
             argc, argv, // Appends --config etc.
             manual_add_config_arguments<ManualConfigArguments{
                 .has_config = true, .requires_input_file = true, .requires_output_file = true}>,
             read_arguments_from_config_file<CreateWBIMDatasetParams>)
      .and_then([](CreateWBIMDatasetParams params) -> rfl::Result<CreateWBIMDatasetParams> {
        if (*params.common->lambda_in <= 0 && *params.lambda_boost <= 0) {
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

auto WIMSketchingExperimentParams::parse_from_args(int argc, char** argv) noexcept
    -> rfl::Result<WIMSketchingExperimentParams> {
  return parse_from_args_generic<WIMSketchingExperimentParams>(
             argc, argv, // Appends --config etc.
             manual_add_config_arguments<ManualConfigArguments{.has_config = true, .requires_input_file = true}>,
             read_arguments_from_config_file<WIMSketchingExperimentParams>)
      .and_then([](WIMSketchingExperimentParams params) -> rfl::Result<WIMSketchingExperimentParams> {
        auto err_msg = sort_and_check_sketching_params(*params.sketching);
        if (err_msg) {
          return rfl::Error{std::move(*err_msg)};
        }
        return std::move(params);
      });
}

auto WBIMSketchingExperimentParams::parse_from_args(int argc, char** argv) noexcept
    -> rfl::Result<WBIMSketchingExperimentParams> {
  return parse_from_args_generic<WBIMSketchingExperimentParams>(
             argc, argv, // Appends --config etc.
             manual_add_config_arguments<ManualConfigArguments{.has_config = true, .requires_input_file = true}>,
             read_arguments_from_config_file<WBIMSketchingExperimentParams>)
      .and_then([](WBIMSketchingExperimentParams params) -> rfl::Result<WBIMSketchingExperimentParams> {
        auto err_msg = sort_and_check_sketching_params(*params.sketching);
        if (err_msg) {
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
        auto err_msg = sort_and_check_sketching_params(*params.sketching);
        if (err_msg) {
          return rfl::Error{std::move(*err_msg)};
        }
        return std::move(params);
      });
}

auto WBIMCoarseningExperimentParams::parse_from_args(int argc, char** argv) noexcept
    -> rfl::Result<WBIMCoarseningExperimentParams> {
  return parse_from_args_generic<WBIMCoarseningExperimentParams>(
             argc, argv, // Appends --config etc.
             manual_add_config_arguments<ManualConfigArguments{.has_config = true, .requires_input_file = true}>,
             read_arguments_from_config_file<WBIMCoarseningExperimentParams>)
      .and_then([](WBIMCoarseningExperimentParams params) -> rfl::Result<WBIMCoarseningExperimentParams> {
        auto err_msg = sort_and_check_sketching_params(*params.sketching);
        if (err_msg) {
          return rfl::Error{std::move(*err_msg)};
        }
        return std::move(params);
      });
}

auto WIMContrastExperimentParams::parse_from_args(int argc, char** argv) noexcept
    -> rfl::Result<WIMContrastExperimentParams> {
  return parse_from_args_generic<WIMContrastExperimentParams>(
             argc, argv, // Appends --config etc.
             manual_add_config_arguments<ManualConfigArguments{.has_config = true, .requires_input_file = true}>,
             read_arguments_from_config_file<WIMContrastExperimentParams>)
      .and_then([](WIMContrastExperimentParams params) -> rfl::Result<WIMContrastExperimentParams> {
        constexpr auto requirements = CheckListArguments{
            .non_empty_required = true,
            .uniqueness_required = true,
            .positive_values_required = true,
        };
        auto n_seeds = std::span{params.n_seeds};
        if (auto err_msg = sort_and_check_list<requirements>("# of seed vertices", n_seeds); err_msg) {
          return rfl::Error{std::move(*err_msg)};
        }
        if (params.with_rr_sketch) {
          auto n_sketches = std::span{params.n_sketches};
          if (auto err_msg = sort_and_check_list<requirements>("# of RR-sketches", n_sketches); err_msg) {
            return rfl::Error{std::move(*err_msg)};
          }
        }
        return std::move(params);
      });
}
