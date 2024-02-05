#include "create_dataset.h"
#include "to_matrix_market.h"
#include "utils/argparse_helper.h"
#include "utils/reflect_for_each_field.h"
#include "wim.h"
#include <fstream>

auto ReadGraphParams::parse_from_args(int argc, char** argv) noexcept -> rfl::Result<ReadGraphParams> {
  return parse_from_args_generic<ReadGraphParams>(
             argc, argv,
             [](ArgumentParser& parser) {
               parser.add_argument("-c", "--config").help("Configuration file as JSON format");
               parser.at("--input-file").required();
               parser.at("--output-file").required();
             },
             [](ArgumentParser& parser, ReadGraphParams& value) {
               if (auto in_path = parser.present("--config")) {
                 read_fields_from_json(value, json::parse(std::ifstream{*in_path}));
               }
             })
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

auto WIMParams::parse_from_args(int argc, char** argv) noexcept -> rfl::Result<WIMParams> {
  return parse_from_args_generic<WIMParams>(
             argc, argv,
             [](ArgumentParser& parser) {
               parser.add_argument("-c", "--config").help("Configuration file as JSON format");
               parser.at("--input-file").required();
             },
             [](ArgumentParser& parser, WIMParams& value) {
               if (auto in_path = parser.present("--config")) {
                 read_fields_from_json(value, json::parse(std::ifstream{*in_path}));
               }
             })
      .and_then([](WIMParams params) -> rfl::Result<WIMParams> {
        auto do_sort_and_check = [](std::string_view name, auto& values) -> ResultVoid {
          if (values.empty()) {
            return rfl::Error{fmt::format("{} can not be an empty list.", name)};
          }
          ranges::sort(values);
          if (auto n0 = values.front(); n0 <= 0) {
            return rfl::Error{fmt::format("{} must be positive integers, while {} is given.", name, n0)};
          }
          constexpr auto adjacent_equal_filter = views::filter(LAMBDA_1(get<0>(_1) == get<1>(_1)));
          for (auto [n0, n1] : views::adjacent<2>(values) | adjacent_equal_filter) {
            return rfl::Error{fmt::format("Duplicated {} (detected {}) is disallowed.", name, n0)};
          }
          return RESULT_VOID_SUCCESS;
        };
        return do_sort_and_check("# of RR-sketches", params.num_sketches)
            .and_then([&](auto) { return do_sort_and_check("# of seed vertices", params.num_seeds); })
            .and_then([&](auto) { return rfl::Result{std::move(params)}; });
      });
}

auto ToMatrixMarketParams::parse_from_args(int argc, char** argv) noexcept -> rfl::Result<ToMatrixMarketParams> {
  return parse_from_args_generic<ToMatrixMarketParams>(
             argc, argv,
             [](ArgumentParser& parser) {
               parser.add_argument("-c", "--config").help("Configuration file as JSON format");
               parser.at("--input-file").required();
               parser.at("--output-file").required();
             },
             [](ArgumentParser& parser, ToMatrixMarketParams& value) {
               if (auto in_path = parser.present("--config")) {
                 read_fields_from_json(value, json::parse(std::ifstream{*in_path}));
               }
             })
      .and_then([](ToMatrixMarketParams params) -> rfl::Result<ToMatrixMarketParams> {
        if (params.is_directed && params.is_undirected) {
          return rfl::Error{"--is-directed and --is-undirected shall not be set at the same time."};
        }
        return std::move(params);
      });
}
