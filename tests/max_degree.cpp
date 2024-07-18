#define BOOST_TEST_MODULE "Max-Degree IM Algorithm"
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
  // Max-degree
  const auto max_degree_ans = std::vector<vertex_id_t>{4, 9, 0, 1, 2, 3, 5, 7, 6, 8};
  BOOST_CHECK_EQUAL_RANGE(max_out_degree(graph), max_degree_ans);
  BOOST_CHECK_EQUAL_RANGE(max_out_degree(inv_graph), max_degree_ans);
  // Max-strength
  const auto max_strength_ans = std::vector<vertex_id_t>{4, 9, 3, 1, 0, 5, 7, 2, 6, 8};
  BOOST_CHECK_EQUAL_RANGE(max_out_strength(graph), max_strength_ans);
  BOOST_CHECK_EQUAL_RANGE(max_out_strength(inv_graph), max_strength_ans);
}

int main(int argc, char* argv[], char* envp[]) {
  easylog::set_min_severity(easylog::Severity::TRACE);
  return boost::unit_test::unit_test_main(init_unit_test, argc, argv);
}
