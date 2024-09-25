#define BOOST_TEST_MODULE Graph Connectivity

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
  auto [n_sccs, component_index] = strongly_connected_components(adj_list);
  std::cout << dump_array("SCC component_index[] of graph1", component_index) << std::endl;
  BOOST_CHECK_EQUAL(n_sccs, 4);
}

BOOST_AUTO_TEST_CASE(graph2) {
  const auto EDGE_LIST_2 = std::vector<EdgeTuple>{{0, 1}, {1, 2}, {2, 3}, {3, 1}};
  auto [adj_list, inv_adj_list] = make_sample_graph(std::span{EDGE_LIST_2});
  // WCC
  BOOST_CHECK_EQUAL(n_weakly_connected_components(adj_list, inv_adj_list), 1);
  BOOST_CHECK_EQUAL(n_weakly_connected_components(adj_list, inv_adj_list, {3}), 1);
  // SCC
  auto [n_sccs, component_index] = strongly_connected_components(adj_list);
  std::cout << dump_array("SCC component_index[] of graph2", component_index) << std::endl;
  BOOST_CHECK_EQUAL(n_sccs, 2);
  auto filter_fn_2_1 = [](vertex_id_t u, vertex_id_t v) {
    return u < v; // Filters out the directed edge (3, 1)
  };
  BOOST_CHECK_EQUAL(n_strongly_connected_components(adj_list, filter_fn_2_1), 4);
}

BOOST_AUTO_TEST_CASE(graph3) {
  const auto EDGE_LIST_3 = std::vector<EdgeTuple>{{0, 1}, {0, 2}, {1, 2}, {1, 4}, {3, 1}, {3, 2}, {4, 3}};
  auto [adj_list, inv_adj_list] = make_sample_graph(std::span{EDGE_LIST_3});
  // WCC
  BOOST_CHECK_EQUAL(n_weakly_connected_components(adj_list, inv_adj_list), 1);
  BOOST_CHECK_EQUAL(n_weakly_connected_components(adj_list, inv_adj_list, {1}), 1);
  BOOST_CHECK_EQUAL(n_weakly_connected_components(adj_list, inv_adj_list, {1, 3}), 2);
  // SCC
  auto [n_sccs, component_index] = strongly_connected_components(adj_list);
  std::cout << dump_array("SCC component_index[] of graph3", component_index) << std::endl;
  BOOST_CHECK_EQUAL(n_sccs, 3);
  auto filter_fn_3_1 = [](vertex_id_t u, vertex_id_t v) {
    return u != 3; // Filters out (3, 1) and (3, 2)
  };
  BOOST_CHECK_EQUAL(n_strongly_connected_components(adj_list, filter_fn_3_1), 5);
}

BOOST_AUTO_TEST_CASE(graph4) {
  const auto EDGE_LIST_4 = std::vector<EdgeTuple>{{0, 1}, {0, 2}, {1, 2}, {1, 4}, {2, 3}, {3, 0}, {3, 1}, {4, 1}};
  auto [adj_list, inv_adj_list] = make_sample_graph(std::span{EDGE_LIST_4});
  // WCC
  BOOST_CHECK_EQUAL(n_weakly_connected_components(adj_list, inv_adj_list), 1);
  BOOST_CHECK_EQUAL(n_weakly_connected_components(adj_list, inv_adj_list, {1}), 2);
  // SCC
  auto [n_sccs, component_index] = strongly_connected_components(adj_list);
  std::cout << dump_array("SCC component_index[] of graph4", component_index) << std::endl;
  BOOST_CHECK_EQUAL(n_sccs, 1);
  auto filter_fn_4_1 = [](vertex_id_t u, vertex_id_t v) {
    return u != 3 || v != 0; // Filters out (3, 0)
  };
  BOOST_CHECK_EQUAL(n_strongly_connected_components(adj_list, filter_fn_4_1), 2);
}

BOOST_AUTO_TEST_CASE(graph5) {
  const auto EDGE_LIST_5 = std::vector<EdgeTuple>{{1, 2}, {1, 3}, {1, 5}, {2, 4}, {3, 1}, {3, 5}, {4, 2}, {5, 1}};
  auto [adj_list, inv_adj_list] = make_sample_graph(std::span{EDGE_LIST_5});
  // WCC
  BOOST_CHECK_EQUAL(n_weakly_connected_components(adj_list, inv_adj_list), 2);
  BOOST_CHECK_EQUAL(n_weakly_connected_components(adj_list, inv_adj_list, {2}), 3);
  BOOST_CHECK_EQUAL(n_weakly_connected_components(adj_list, inv_adj_list, {1, 5}), 3);
  // SCC
  auto [n_sccs, component_index] = strongly_connected_components(adj_list);
  std::cout << dump_array("SCC component_index[] of graph5", component_index) << std::endl;
  BOOST_CHECK_EQUAL(n_sccs, 3);
  auto filter_fn_5_1 = [](vertex_id_t u, vertex_id_t v) {
    std::tie(u, v) = std::minmax(u, v);
    return u != 1 || v != 3; // Filters out (1, 3) and (3, 1)
  };
  BOOST_CHECK_EQUAL(n_strongly_connected_components(adj_list, filter_fn_5_1), 4);
}
