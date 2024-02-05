#pragma once

#include "utils/result.h"

struct ToMatrixMarketParams {
  std::string input_file;
  std::string output_file;

  unsigned num_lines_ignored = 0;
  bool is_directed = false;
  bool is_undirected = false;

  static auto parse_from_args(int argc, char** argv) noexcept -> rfl::Result<ToMatrixMarketParams>;
};

auto graph_text_to_matrix_market(const ToMatrixMarketParams& params) -> ResultVoid;
