#include "to_matrix_market.h"
#include "utils/boost_assert.h"
#include "utils/graph.h"
#include "utils/utils.h"
#include <fmt/format.h>
#include <fstream>
#include <map>
#include <ylt/easylog.hpp>

namespace {
constexpr auto NUM_PREVIEW_LINES = 10;

auto dump_as_matrix_market_text(size_t n, std::span<const std::pair<vertex_id_t, vertex_id_t>> edges,
                                bool is_undirected) {
  auto m = ranges::size(edges);
  auto res = fmt::format("%%MatrixMarket matrix coordinate pattern {0}\n{1} {1} {2}\n",
                         (is_undirected ? "symmetric" : "general"), n, m);
  // Shifts 0-based to 1-based in MatrixMarket
  for (auto [u, v] : edges) {
    res += fmt::format("{} {}\n", u + 1, v + 1);
  }
  return res;
}

auto first_lines(std::string_view str, size_t n, char delim = '\n') -> std::string_view {
  auto pos = 0zu;
  for (auto i : range(n)) {
    pos = str.find(delim, pos);
    if (pos == std::string_view::npos) {
      return str;
    }
    pos += 1;
  }
  return str.substr(0, pos - 1); // -1 : Excludes the trailing '\n'
}
} // namespace

auto graph_text_to_matrix_market(const ToMatrixMarketParams& params) -> ResultVoid {
  auto fin = std::ifstream{params.input_file};
  if (!fin.is_open()) {
    return rfl::Error{fmt::format("Failed to open input file '{}'.", params.input_file)};
  }
  auto fout = std::ofstream{params.output_file};
  if (!fout.is_open()) {
    return rfl::Error{fmt::format("Failed to open output file '{}'.", params.output_file)};
  }

  auto index_map = std::map<vertex_id_t, vertex_id_t>{};
  auto mapped_edges = std::vector<std::pair<vertex_id_t, vertex_id_t>>{};
  auto self_loop_count = uint64_t{0};
  auto ascending_flag = false;
  auto descending_flag = false;

  auto get_mapped_index = [&](vertex_id_t v) {
    auto it = index_map.find(v);
    if (it != index_map.end()) {
      return it->second;
    }
    auto res = static_cast<vertex_id_t>(index_map.size());
    index_map.emplace(v, res);
    return res;
  };

  // Ignores the first several lines
  if (params.num_lines_ignored > 0) {
    constexpr auto msg_pattern = "Ignores the first {} line(s) of input file '{}'.";
    ELOGFMT(INFO, msg_pattern, params.num_lines_ignored, params.input_file);
    for (auto i : range(params.num_lines_ignored)) {
      fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
  }
  for (vertex_id_t u, v; fin >> u >> v;) {
    if (u == v) {
      self_loop_count += 1;
      continue; // Skips self-loops
    } else if (u < v) {
      ascending_flag = true;
    } else {
      descending_flag = true;
    }
    mapped_edges.emplace_back(get_mapped_index(u), get_mapped_index(v));
  }
  if (self_loop_count > 0) {
    ELOGFMT(WARNING, "{} self-loops detected & removed.", self_loop_count);
  }
  // Undirected if all the conditions below are satisfied:
  // (1) --is-undirected is set;
  // (2) only one of ascending_flag (i.e. each edge (u, v) satisfied u < v)
  //     or descending_flag (i.e. each edge (u, v) satisfied u > v) is enabled.
  auto is_undirected = [&] {
    BOOST_ASSERT_MSG(!params.is_directed || !params.is_undirected, "Two flags can't be set at the same time.");
    if (params.is_undirected) {
      ELOG_INFO << "--is-undirected flag is set.";
      return true;
    } else if (params.is_directed) {
      ELOG_INFO << "--is-directed flag is set.";
      return false;
    } else if (!ascending_flag || !descending_flag) {
      ELOG_INFO << "Detected as undirected graph since every edge (u, v) satisfies "
                << (ascending_flag ? "u < v" : "u > v") << ".";
      return true;
    } else {
      ELOG_INFO << "Detected as directed graph.";
      return false;
    }
  }();
  // Sorts & removes duplicated edges
  ranges::sort(mapped_edges);
  auto [erase_begin, erase_end] = ranges::unique(mapped_edges);
  if (erase_begin != erase_end) {
    ELOGFMT(WARNING, "{} duplicated edges detected & removed.", erase_end - erase_begin);
    mapped_edges.erase(erase_begin, erase_end);
  }
  // Output the MatrixMarket text file
  auto mm_text = dump_as_matrix_market_text(index_map.size(), mapped_edges, is_undirected);
  ELOGFMT(DEBUG, "The first {} lines of Matrix Market output:\n{}", NUM_PREVIEW_LINES,
          first_lines(mm_text, NUM_PREVIEW_LINES));
  fout.write(mm_text.c_str(), mm_text.length());

  return RESULT_VOID_SUCCESS;
}
