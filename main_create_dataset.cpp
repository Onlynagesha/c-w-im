#include "create_dataset.h"
#include "dump.h"
#include <fmt/ranges.h>

#define VERIFY_CHECK(expr1, expr2)                                      \
  do {                                                                  \
    auto check_result = (expr1) == (expr2);                             \
    if (check_result) {                                                 \
      fmt::println("PASS: {} == {}", #expr1, #expr2);                   \
    } else {                                                            \
      return rfl::Error{fmt::format("FAIL: {} != {}", #expr1, #expr2)}; \
    }                                                                   \
  } while (false)

constexpr auto VERIFY_NUM_EDGES_THRESHOLD = 20;
constexpr auto PREVIEW_COUNT = 10;

auto verify_correctness(const std::string& data_file, const DirectedEdgeList<WIMEdge>& expected_graph) -> ResultVoid {
  return read_directed_wim_edge_list(data_file).and_then([&](const auto& read_result) -> ResultVoid {
    fmt::println("Starts verification.");
    // Checking |V| and |E|
    VERIFY_CHECK(read_result.edge_list.num_vertices(), expected_graph.num_vertices());
    VERIFY_CHECK(read_result.edge_list.num_edges(), expected_graph.num_edges());
    // Checking the consistency of each edge
    for (auto [i, e1, e2] : views::zip(views::iota(0), read_result.edge_list, expected_graph)) {
      if (e1 != e2) {
        return rfl::Error{fmt::format("FAIL since edge #{} differs: {} vs. {}", i, e1, e2)};
      }
    }
    fmt::println("PASS: Equality check of each edge.");
    return RESULT_VOID_SUCCESS;
  });
}

int main(int argc, char** argv) {
  auto params = ReadGraphParams::parse_from_args(argc, argv);
  if (!params) {
    fmt::println(stderr, "Failed to parse parameters: `{}'", params.error()->what());
    return -1;
  }
  auto graph = create_wim_dataset(*params);
  fmt::println("Preview of first (up to) {} vertices:", PREVIEW_COUNT);
  for (auto u : vertices(graph.edge_list) | views::take(PREVIEW_COUNT)) {
    fmt::println("\tc({}) = {}", u, graph.vertex_weights[u]);
  }
  fmt::println("Preview of first (up to) {} edges:", PREVIEW_COUNT);
  for (auto [u, v, w] : graph.edge_list | views::take(PREVIEW_COUNT)) {
    fmt::println("\t{} -> {}: properties = {}", u, v, w);
  }

  auto res = write_directed_wim_edge_list_r(graph, (*params).output_file).and_then([&](auto) -> ResultVoid {
    fmt::println("Successfully serializes dataset from input '{}' to output '{}'.", (*params).input_file,
                 (*params).output_file);
    if (graph.edge_list.num_edges() <= VERIFY_NUM_EDGES_THRESHOLD) {
      return verify_correctness((*params).output_file, graph.edge_list);
    }
    return RESULT_VOID_SUCCESS;
  });
  if (!res) {
    constexpr auto msg_pattern = "Failed to write the dataset to file '{}': `{}'";
    fmt::println(stderr, msg_pattern, (*params).output_file, res.error()->what());
    return -1;
  }
  return 0;
}
