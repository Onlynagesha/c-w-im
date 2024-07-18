#define BOOST_TEST_MODULE "IMRank Algorithm"
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
  // Expected Mr(0): [1.1188, 0.8426, 1.1232, 0.4639, 2.1512, 0.5576, 0.4480, 0.6540, 0.3510, 2.2899]
  const auto ans = std::vector<vertex_id_t>{9, 4, 2, 0, 7, 1, 5, 3, 6, 8};
  auto out = wim_imrank(inv_graph, {.n_iterations = 10, .n_iterations_before_topk_fixed = 3});
  BOOST_CHECK_EQUAL_RANGE(out, ans);
}

int main(int argc, char* argv[], char* envp[]) {
  easylog::set_min_severity(easylog::Severity::TRACE);
  return boost::unit_test::unit_test_main(init_unit_test, argc, argv);
}
