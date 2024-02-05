#include "create_dataset.h"
#include "dump.h"
#include <fmt/ranges.h>
#include <ylt/easylog.hpp>

#define VERIFY_CHECK(expr1, expr2)                                      \
  do {                                                                  \
    auto check_result = (expr1) == (expr2);                             \
    if (check_result) {                                                 \
      ELOGFMT(INFO, "PASS: {} == {}", #expr1, #expr2);                  \
    } else {                                                            \
      return rfl::Error{fmt::format("FAIL: {} != {}", #expr1, #expr2)}; \
    }                                                                   \
  } while (false)

constexpr auto VERIFY_NUM_EDGES_THRESHOLD = 20;
constexpr auto PREVIEW_COUNT = 10;

auto verify_correctness(const std::string& data_file, const DirectedEdgeList<WIMEdge>& expected_graph) -> ResultVoid {
  return read_directed_wim_edge_list(data_file).and_then([&](const auto& read_result) -> ResultVoid {
    ELOG_INFO << "Starts verification.";
    // Checking |V| and |E|
    VERIFY_CHECK(read_result.edge_list.num_vertices(), expected_graph.num_vertices());
    VERIFY_CHECK(read_result.edge_list.num_edges(), expected_graph.num_edges());
    // Checking the consistency of each edge
    for (auto [i, e1, e2] : views::zip(views::iota(0), read_result.edge_list, expected_graph)) {
      if (e1 != e2) {
        return rfl::Error{fmt::format("FAIL since edge #{} differs: {} vs. {}", i, e1, e2)};
      }
    }
    ELOG_INFO << "PASS: Equality check of each edge.";
    return RESULT_VOID_SUCCESS;
  });
}

auto main_worker(int argc, char** argv) -> ResultVoid try {
  return ReadGraphParams::parse_from_args(argc, argv).and_then([](ReadGraphParams params) {
    easylog::set_min_severity(params.log_level);
    auto graph = create_wim_dataset(params);

    auto preview_str = fmt::format("Preview of first (up to) {} vertices:", PREVIEW_COUNT);
    for (auto u : vertices(graph.edge_list) | views::take(PREVIEW_COUNT)) {
      preview_str += fmt::format("\n\tc({}) = {:.4}", u, graph.vertex_weights[u]);
    }
    preview_str += fmt::format("\nPreview of first (up to) {} edges:", PREVIEW_COUNT);
    for (auto [u, v, w] : graph.edge_list | views::take(PREVIEW_COUNT)) {
      preview_str += fmt::format("\n\t{} -> {}: properties = {}", u, v, w);
    }
    ELOG_DEBUG << preview_str;

    return write_directed_wim_edge_list_r(graph, params.output_file).and_then([&](auto) -> ResultVoid {
      constexpr auto msg_pattern = "Successfully serializes dataset from input '{}' to output '{}'.";
      ELOGFMT(INFO, msg_pattern, params.input_file, params.output_file);
      // Verification for small graphs
      return (graph.edge_list.num_edges() <= VERIFY_NUM_EDGES_THRESHOLD)
                 ? verify_correctness(params.output_file, graph.edge_list)
                 : RESULT_VOID_SUCCESS;
    });
  });
}
RFL_RESULT_CATCH_HANDLER()

int main(int argc, char** argv) {
  return main_worker(argc, argv)
      .and_then([](auto) {
        ELOG_INFO << "Done creating dataset.";
        return rfl::Result{0};
      })
      .or_else([](const rfl::Error& error) {
        ELOGFMT(CRITICAL, "ERROR: `{}'", error.what());
        return rfl::Result{-1};
      })
      .value_or(-1);
}
