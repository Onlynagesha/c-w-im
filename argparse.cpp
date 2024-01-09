#include "create_dataset.h"
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
      });
}

auto WIMParams::parse_from_args(int argc, char** argv) noexcept -> rfl::Result<WIMParams> {
  auto res = parse_from_args_generic<WIMParams>(
      argc, argv,
      [](ArgumentParser& parser) {
        parser.add_argument("-c", "--config").help("Configuration file as JSON format");
        parser.at("--input-file").required();
      },
      [](ArgumentParser& parser, WIMParams& value) {
        if (auto in_path = parser.present("--config")) {
          read_fields_from_json(value, json::parse(std::ifstream{*in_path}));
        }
      });
  return res.and_then([](WIMParams params) -> rfl::Result<WIMParams> {
    // (1) Checks & sorts num_sketches
    if (params.num_sketches.empty()) {
      return rfl::Error{"# of RR-sketches can not be an empty list."};
    }
    ranges::sort(params.num_sketches);
    if (auto n0 = params.num_sketches.front(); n0 <= 0) {
      return rfl::Error{fmt::format("# of RR-sketches must be positive integers, while {} is given.", n0)};
    }
    auto adjacent_equal_filter = views::filter(LAMBDA_1(get<0>(_1) == get<1>(_1)));
    for (auto [n0, _] : views::adjacent<2>(params.num_sketches) | adjacent_equal_filter) {
      return rfl::Error{fmt::format("Duplicated # of RR-sketches (detected {}) is disallowed.", n0)};
    }
    // (2) Checks & sorts num_seeds
    if (params.num_seeds.empty()) {
      return rfl::Error{"# of seed vertices can not be an empty list."};
    }
    ranges::sort(params.num_seeds);
    if (auto n0 = params.num_seeds.front(); n0 <= 0) {
      return rfl::Error{fmt::format("# of seed vertices must be positive integers, while {} is given.", n0)};
    }
    for (auto [n0, _] : views::adjacent<2>(params.num_seeds) | adjacent_equal_filter) {
      return rfl::Error{fmt::format("Duplicated # of seeds (detected {}) is disallowed.", n0)};
    }
    // Done all checking
    return std::move(params);
  });
}
