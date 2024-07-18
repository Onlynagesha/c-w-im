#define BOOST_TEST_MODULE "Pagerank IM Algorithm"
#define BOOST_TEST_DYN_LINK

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API

#include "contrast_algorithms.h"
#include "tests/sample_graph.h"
#include "utils/easylog.h"
#include <boost/test/unit_test.hpp>
#include <fmt/ranges.h>

#define BOOST_CHECK_EQUAL_RANGE(lhs, rhs) BOOST_CHECK_EQUAL(fmt::format("{}", lhs), fmt::format("{}", rhs))

BOOST_AUTO_TEST_CASE(graphA) {
  auto [graph, inv_graph] = make_sample_wim_graph("A");
  auto [n, m] = graph_n_m(graph);
  auto vertex_weights = [&]() {
    auto view = views::iota(vertex_id_t{10}, vertex_id_t{10} + graph::num_vertices(graph));
    return std::vector<vertex_weight_t>(view.begin(), view.end());
  }();

  struct TestParams {
    bool transpose = false;
    bool weighted = false;
  };
  auto do_max_pagerank = [&](TestParams test_params = TestParams{}) {
    auto algo_params = PagerankParams{
        .k = n,
        .n_iterations = 40,
        .transpose = test_params.transpose,
        .uses_vertex_weight = test_params.weighted,
        .uses_edge_weight = test_params.weighted,
    };
    return max_pagerank(graph, inv_graph, vertex_weights, algo_params);
  };

  using AnsILWrapper = std::initializer_list<std::vector<vertex_id_t>>;

  // Expected PR(1): [0.8583, 0.7875, 0.5750, 1.0708, 1.8500, 0.7875, 1.7083, 1.2125, 1.0000, 0.1500]
  const auto ans = std::vector<vertex_id_t>{7, 4, 6, 8, 3, 5, 1, 2, 0, 9};
  BOOST_CHECK_EQUAL_RANGE(do_max_pagerank(), AnsILWrapper{ans});

  // Expected PR(1): [10.6023, 7.2249, 5.7705, 14.4235, 21.8008, 13.7980, 27.6863, 20.6125, 20.2313, 2.8500]
  const auto ans_weighted = std::vector<vertex_id_t>{7, 6, 8, 4, 5, 3, 1, 2, 0, 9};
  BOOST_CHECK_EQUAL_RANGE(do_max_pagerank({.weighted = true}), AnsILWrapper{ans_weighted});

  // Expected PR(1): [0.8583, 1.2125, 1.0000, 0.6458, 1.7083, 0.7875, 0.5750, 0.7875, 0.4333, 1.1417]
  const auto ans_transpose = std::vector<vertex_id_t>{4, 1, 0, 9, 2, 5, 7, 3, 6, 8};
  BOOST_CHECK_EQUAL_RANGE(do_max_pagerank({.transpose = true}), AnsILWrapper{ans_transpose});

  // Expected PR(1): [11.9815, 16.0856, 10.7138, 13.1558, 22.2212, 13.2584, 8.6190, 11.2789, 4.6947, 16.8412]
  const auto ans_transpose_weighted = std::vector<vertex_id_t>{9, 4, 1, 0, 3, 2, 5, 7, 6, 8};
  BOOST_CHECK_EQUAL_RANGE(do_max_pagerank({.transpose = true, .weighted = true}), AnsILWrapper{ans_transpose_weighted});
}

int main(int argc, char* argv[], char* envp[]) {
  easylog::set_min_severity(easylog::Severity::TRACE);
  return boost::unit_test::unit_test_main(init_unit_test, argc, argv);
}
