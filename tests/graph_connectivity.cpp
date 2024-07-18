#define BOOST_TEST_MODULE Graph Connectivity
#define BOOST_TEST_DYN_LINK

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API

#include "graph_connectivity.h"
#include "tests/sample_graph.h"
#include "utils/easylog.h"
#include <boost/test/unit_test.hpp>

using EdgeTuple = std::tuple<vertex_id_t, vertex_id_t>;

BOOST_AUTO_TEST_CASE(graph1) {
  const auto EDGE_LIST_1 = std::vector<EdgeTuple>{{0, 1}, {1, 2}, {1, 3}, {2, 3}};
  auto [adj_list, inv_adj_list] = make_sample_graph(std::span{EDGE_LIST_1});
  // WCC
  BOOST_CHECK_EQUAL(n_weakly_connected_components(adj_list, inv_adj_list), 1);
  BOOST_CHECK_EQUAL(n_weakly_connected_components(adj_list, inv_adj_list, {1}), 2);
  // SCC
  BOOST_CHECK_EQUAL(n_strongly_connected_components(adj_list), 4);
}

BOOST_AUTO_TEST_CASE(graph2) {
  const auto EDGE_LIST_2 = std::vector<EdgeTuple>{{0, 1}, {1, 2}, {2, 3}, {3, 1}};
  auto [adj_list, inv_adj_list] = make_sample_graph(std::span{EDGE_LIST_2});
  // WCC
  BOOST_CHECK_EQUAL(n_weakly_connected_components(adj_list, inv_adj_list), 1);
  BOOST_CHECK_EQUAL(n_weakly_connected_components(adj_list, inv_adj_list, {3}), 1);
  // SCC
  BOOST_CHECK_EQUAL(n_strongly_connected_components(adj_list), 2);
}

BOOST_AUTO_TEST_CASE(graph3) {
  const auto EDGE_LIST_3 = std::vector<EdgeTuple>{{0, 1}, {0, 2}, {1, 2}, {1, 4}, {3, 1}, {3, 2}, {4, 3}};
  auto [adj_list, inv_adj_list] = make_sample_graph(std::span{EDGE_LIST_3});
  // WCC
  BOOST_CHECK_EQUAL(n_weakly_connected_components(adj_list, inv_adj_list), 1);
  BOOST_CHECK_EQUAL(n_weakly_connected_components(adj_list, inv_adj_list, {1}), 1);
  BOOST_CHECK_EQUAL(n_weakly_connected_components(adj_list, inv_adj_list, {1, 3}), 2);
  // SCC
  BOOST_CHECK_EQUAL(n_strongly_connected_components(adj_list), 3);
}

BOOST_AUTO_TEST_CASE(graph4) {
  const auto EDGE_LIST_4 = std::vector<EdgeTuple>{{0, 1}, {0, 2}, {1, 2}, {1, 4}, {2, 3}, {3, 0}, {3, 1}, {4, 1}};
  auto [adj_list, inv_adj_list] = make_sample_graph(std::span{EDGE_LIST_4});
  // WCC
  BOOST_CHECK_EQUAL(n_weakly_connected_components(adj_list, inv_adj_list), 1);
  BOOST_CHECK_EQUAL(n_weakly_connected_components(adj_list, inv_adj_list, {1}), 2);
  // SCC
  BOOST_CHECK_EQUAL(n_strongly_connected_components(adj_list), 1);
}

BOOST_AUTO_TEST_CASE(graph5) {
  const auto EDGE_LIST_5 = std::vector<EdgeTuple>{{1, 2}, {1, 3}, {1, 5}, {2, 4}, {3, 1}, {3, 5}, {4, 2}, {5, 1}};
  auto [adj_list, inv_adj_list] = make_sample_graph(std::span{EDGE_LIST_5});
  // WCC
  BOOST_CHECK_EQUAL(n_weakly_connected_components(adj_list, inv_adj_list), 2);
  BOOST_CHECK_EQUAL(n_weakly_connected_components(adj_list, inv_adj_list, {2}), 3);
  BOOST_CHECK_EQUAL(n_weakly_connected_components(adj_list, inv_adj_list, {1, 5}), 3);
  // SCC
  BOOST_CHECK_EQUAL(n_strongly_connected_components(adj_list), 3);
}

int main(int argc, char* argv[], char* envp[]) {
  easylog::set_min_severity(easylog::Severity::TRACE);
  return boost::unit_test::unit_test_main(init_unit_test, argc, argv);
}
