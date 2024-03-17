#include "coarsening.h"
#include "dump.h"
#include "playground/sample_graph.h"
#include "utils/easylog.h"
#include "wim.h"
#include <fmt/ranges.h>
#include <nwgraph/adaptors/edge_range.hpp>

int main() {
  easylog::set_min_severity(easylog::Severity::TRACE);
  easylog::set_async(false);

  auto [graph, inv_graph] = make_sample_wbim_graph_1();
  auto [n, m] = graph_n_m(graph);
  auto vertex_weights = [&]() {
    auto view = views::iota(vertex_id_t{10}, vertex_id_t{10} + n);
    return std::vector<vertex_weight_t>(view.begin(), view.end());
  }();

  auto coarsening_params = CoarseningParams{
      .neighbor_match_rule = NeighborMatchRule::HEM_P_MAX,
      .edge_weight_rule = EdgeWeightRule::SEPARATE_SIMPLE,
      .boosted_edge_weight_rule = BoostedEdgeWeightRule::BEST_BOOSTED_INDEX,
      .boosted_selection_rule = BoostedSelectionRule::AS_SOURCE,
      .in_out_heuristic_rule = InOutHeuristicRule::P,
      .vertex_weight_rule = VertexWeightRule::AVERAGE_BY_PATHS,
      .max_distance_from_seed = 10,
  };
  auto expanding_params = ExpandingParams{
      .vertex_expanding_rule = VertexExpandingRule::ITERATIVE, .n_iterations = 5, .simulation_try_count = 5};

  auto bidir_graph = merge_wbim_edge_to_undirected(graph, coarsening_params);
  ELOG_INFO << [&] {
    constexpr auto msg_pattern_header = "Merged bidirectional graph: |V|, |E| = {}";
    auto res = fmt::format(msg_pattern_header, graph_n_m(bidir_graph));
    for (auto [u, v, p] : graph::make_edge_range<0>(bidir_graph)) {
      res += fmt::format("\n\tu = {}, v = {}, p = {}", u, v, p);
    }
    return res;
  }();

  auto n_groups = 5;
  auto group_id = std::vector<vertex_id_t>{0, 4, 2, 0, 1, 2, 0, 1, 2, 3};
  auto seeds = VertexSet(n, {1, 9});
  // auto [n_groups, group_id] = mongoose_match(bidir_graph, coarsening_params);

  auto detailed_res = coarsen_wbim_graph_by_match_d( //
      graph, inv_graph, vertex_weights, seeds, n_groups, group_id, coarsening_params);
  ELOGFMT(INFO, "Detailed coarsening result: {:4}", *detailed_res);

  return 0;
}

/*
"AS_SOURCE": Best boosted vertex selection:
estimated gain for group [0, 3, 6] = [0.0963, 0.2442, 0.1054] Best boosted vertex selection: estimated gain for group
[4, 7] = [0.6863, 0.1583] Best boosted vertex selection: estimated gain for group [2, 5, 8] = [0.0684, 0.0926, 0.0385]

"AS_TARGET":
Best boosted vertex selection: estimated gain for group [0, 3, 6] = [0.0000, 0.2320, 0.2139]
Best boosted vertex selection: estimated gain for group [4, 7] = [0.5728, 0.2718]
Best boosted vertex selection: estimated gain for group [2, 5, 8] = [0.0310, 0.1066, 0.0619]

Starts WBIM coarsening: n = 10, n_groups = 5, parameters: {
    "boosted_edge_weight_rule": "BEST_BOOSTED_INDEX",
    "edge_weight_rule": "SEPARATE_SIMPLE",
    "in_out_heuristic_rule": "P",
    "max_distance_from_seed": 6,
    "neighbor_match_rule": "HEM_P_MAX",
    "seed_edge_weight_rule": "BEST_SEED_INDEX",
    "vertex_weight_rule": "AVERAGE_BY_PATHS"
}
.coarsened.adj_list = {
    (0, 1): {"p":0.29203855991363525,"p_boost":0.42106086015701294},
    (1, 0): {"p":0.2200000286102295,"p_boost":0.4399999976158142},
    (1, 2): {"p":0.4227813482284546,"p_boost":0.627460777759552},
    (2, 0): {"p":0.09320050477981567,"p_boost":0.09320050477981567},
    (2, 1): {"p":0.09890002012252808,"p_boost":0.19780004024505615},
    (3, 0): {"p":0.6639999747276306,"p_boost":0.9520000219345093},
    (4, 1): {"p":0.41999995708465576,"p_boost":0.8399999737739563},
    (4, 2): {"p":0.3100000023841858,"p_boost":0.6200000047683716}
},
.coarsened.inv_adj_list = {
    (0, 1): {"p":0.2200000286102295,"p_boost":0.4399999976158142},
    (0, 2): {"p":0.09320050477981567,"p_boost":0.09320050477981567},
    (0, 3): {"p":0.6639999747276306,"p_boost":0.9520000219345093},
    (1, 0): {"p":0.29203855991363525,"p_boost":0.42106086015701294},
    (1, 2): {"p":0.09890002012252808,"p_boost":0.19780004024505615},
    (1, 4): {"p":0.41999995708465576,"p_boost":0.8399999737739563},
    (2, 1): {"p":0.4227813482284546,"p_boost":0.627460777759552},
    (2, 4): {"p":0.3100000023841858,"p_boost":0.6200000047683716}
},
.coarsened.vertex_weights = {
    16.5917, 19.1710, 18.8621, 19.0000, 11.0000
},

Starts WBIM coarsening: n = 10, n_groups = 5, parameters: {
    "boosted_edge_weight_rule": "BEST_BOOSTED_INDEX",
    "boosted_selection_rule": "AS_SOURCE",
    "edge_weight_rule": "SEPARATE_SIMPLE",
    "in_out_heuristic_rule": "SEED",
    "max_distance_from_seed": 10,
    "neighbor_match_rule": "HEM_P_MAX",
    "seed_edge_weight_rule": "BEST_SEED_INDEX",
    "vertex_weight_rule": "AVERAGE_BY_PATHS"
}
F[10][]       = [0.4186, 0.0000, 0.3100, 0.4941, 0.5796, 0.2991, 0.3908, 0.3590, 0.2459, 0.0000]
F_boost[10][] = [0.8124, 0.0000, 0.6200, 0.8043, 0.9199, 0.5471, 0.6792, 0.6487, 0.4573, 0.0000]
.coarsened.adj_list = {
    (0, 1): {"p":0.321037620306015,"p_boost":0.46373075246810913},
    (1, 0): {"p":0.15351934731006622,"p_boost":0.30703866481781006},
    (1, 2): {"p":0.42004096508026123,"p_boost":0.6089481711387634},
    (2, 0): {"p":0.09269706904888153,"p_boost":0.09269706904888153},
    (2, 1): {"p":0.11631104350090027,"p_boost":0.23262208700180054},
    (3, 0): {"p":0.6639999747276306,"p_boost":0.8079999685287476},
    (4, 1): {"p":0.41999995708465576,"p_boost":0.8399999737739563},
    (4, 2): {"p":0.3100000023841858,"p_boost":0.6200000047683716}
},
.coarsened.inv_adj_list = {
    (0, 1): {"p":0.15351934731006622,"p_boost":0.30703866481781006},
    (0, 2): {"p":0.09269706904888153,"p_boost":0.09269706904888153},
    (0, 3): {"p":0.6639999747276306,"p_boost":0.8079999685287476},
    (1, 0): {"p":0.321037620306015,"p_boost":0.46373075246810913},
    (1, 2): {"p":0.11631104350090027,"p_boost":0.23262208700180054},
    (1, 4): {"p":0.41999995708465576,"p_boost":0.8399999737739563},
    (2, 1): {"p":0.42004096508026123,"p_boost":0.6089481711387634},
    (2, 4): {"p":0.3100000023841858,"p_boost":0.6200000047683716}
},
.coarsened.vertex_weights = {
    16.4080, 18.4946, 19.1306, 19.0000, 11.0000
},
*/