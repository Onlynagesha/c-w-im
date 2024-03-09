#include "contrast_algorithms.h"
#include "playground/sample_graph.h"
#include "utils/easylog.h"
#include <fmt/ranges.h>

auto test_max_degree() -> void {
  auto [graph, inv_graph] = make_sample_wim_graph_1();
  // Expected result: [4, 9, 0, 1, 2, 3, 5, 7, 6, 8]
  ELOGFMT(INFO, "Max Out-egree = {}", max_out_degree(graph));
  ELOGFMT(INFO, "With transpose graph given: Max out-degree = {}", max_out_degree(inv_graph));
  // Expected result: [4, 9, 3, 1, 0, 5, 7, 2, 6, 8]
  ELOGFMT(INFO, "Max Out-strength = {}", max_out_strength(graph));
  ELOGFMT(INFO, "With transpose graph given: Max out-strength = {}", max_out_strength(inv_graph));
}

auto test_pagerank() -> void {
  auto [graph, inv_graph] = make_sample_wim_graph_1();
  auto vertex_weights = [&]() {
    auto view = views::iota(vertex_id_t{10}, vertex_id_t{10} + graph::num_vertices(graph));
    return std::vector<vertex_weight_t>(view.begin(), view.end());
  }();

  auto do_test_pagerank = [&](uint64_t n_iterations, vertex_id_t top_k, bool transpose) {
    // transpose == false:
    //   Expected PR(1): [0.8583, 0.7875, 0.5750, 1.0708, 1.8500, 0.7875, 1.7083, 1.2125, 1.0000, 0.1500]
    // transpose == true:
    //   Expected PR(1): [0.8583, 1.2125, 1.0000, 0.6458, 1.7083, 0.7875, 0.5750, 0.7875, 0.4333, 1.1417]
    auto result = max_pagerank(graph, inv_graph, vertex_weights,
                               {.k = 3,
                                .n_iterations = n_iterations,
                                .transpose = transpose,
                                .uses_vertex_weight = false,
                                .uses_edge_weight = false});
    ELOGFMT(INFO, "Testing Pagerank done: Result = {}\n", result);

    // transpose == false:
    //   Expected PR(1): [10.6023, 7.2249, 5.7705, 14.4235, 21.8008, 13.7980, 27.6863, 20.6125, 20.2313, 2.8500]
    // transpose == true:
    //   Expected PR(1): [11.9815, 16.0856, 10.7138, 13.1558, 22.2212, 13.2584, 8.6190, 11.2789, 4.6947, 16.8412]
    auto result_v_w = max_pagerank(graph, inv_graph, vertex_weights,
                                   {.k = 3,
                                    .n_iterations = n_iterations,
                                    .transpose = transpose,
                                    .uses_vertex_weight = true,
                                    .uses_edge_weight = true});
    ELOGFMT(INFO, "Testing Pagerank(v, w) done: Result = {}\n\n", result_v_w);
  };
  ELOG_INFO << "Testing Pagerank of G(V, E):";
  do_test_pagerank(100, 5, false);
  ELOG_INFO << "Testing Pagerank of the transpose graph G'(V, E'):";
  do_test_pagerank(100, 5, true);
}

auto test_imrank() -> void {
  ELOG_INFO << "Testing IMRank algorithm.";
  auto [graph, inv_graph] = make_sample_wim_graph_1();
  // Expected Mr(0): [1.1188, 0.8426, 1.1232, 0.4639, 2.1512, 0.5576, 0.4480, 0.6540, 0.3510, 2.2899]
  auto imrank_res = wim_imrank(inv_graph, {.n_iterations = 100, .n_iterations_before_topk_fixed = 3});
  ELOGFMT(INFO, "Result of IMRank algorithm = {}", imrank_res);
}

int main() {
  easylog::set_min_severity(easylog::Severity::TRACE);
  test_max_degree();
  test_pagerank();
  test_imrank();
  return 0;
}
