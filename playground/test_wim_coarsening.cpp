#include "coarsening.h"
#include "dump.h"
#include "playground/sample_graph.h"
#include "utils/easylog.h"
#include "wim.h"
#include <fmt/ranges.h>
#include <nwgraph/adaptors/edge_range.hpp>

auto monte_carlo_simulate(const AdjacencyList<WIMEdge>& graph, const VertexSet& seeds, const VertexSet& dest,
                          uint64_t try_count) -> double {
  auto success_count = 0.0;
  auto queue = std::vector<vertex_id_t>();
  auto vis = DynamicBitset{};
  ELOGFMT(INFO, "seeds = {}, {:b}, dest = {}, {:b}", //
          seeds.vertex_list, seeds.mask.to_ulong(), dest.vertex_list, dest.mask.to_ulong());
  for (auto attempt_index : range(try_count)) {
    queue = {seeds.vertex_list.begin(), seeds.vertex_list.end()};
    vis = seeds.mask;

    for (size_t qi = 0; qi < queue.size(); qi++) {
      auto cur = queue[qi];
      if (dest.contains(cur)) {
        // ELOGFMT(INFO, "Attempt #{}: Current queue = {}", attempt_index, queue);
        success_count += 1.0;
        break;
      }
      for (auto [v, w] : graph[cur]) {
        if (!vis.test(v) && rand_bool(w.p)) {
          queue.push_back(v);
          vis.set(v);
        }
      }
    }
  }
  return success_count / try_count;
}

auto monte_carlo_test(const AdjacencyList<WIMEdge>& graph) {
  ELOG_INFO << [&] {
    constexpr auto msg_pattern_header = "Current graph to be tested: |V|, |E| = {}";
    auto res = fmt::format(msg_pattern_header, graph_n_m(graph));
    for (auto [u, v, p] : graph::make_edge_range<0>(as_non_const(graph))) {
      res += fmt::format("\n\tu = {}, v = {}, p = {}", u, v, p);
    }
    return res;
  }();
  constexpr auto TRY_COUNT = 20'000'000uLL;
  for (auto seed : std::vector<vertex_id_t>{0, 2, 4, 6}) {
    auto result = monte_carlo_simulate(graph, {7, {seed}}, {7, {1, 3, 5}}, TRY_COUNT);
    ELOGFMT(INFO, "Simulated result with seed {} = {:.6f}", seed, result);
  }
  for (auto seeds : std::vector<VertexSet>{{7, {0, 2}}, {7, {0, 4}}, {7, {2, 4}}}) {
    auto result = monte_carlo_simulate(graph, seeds, {7, {1, 3, 5}}, TRY_COUNT);
    ELOGFMT(INFO, "Simulated result with seed {} = {:.6f}", seeds.vertex_list, result);
  }
}

int main() {
  easylog::set_min_severity(easylog::Severity::TRACE);
  easylog::set_async(false);

  auto [graph, inv_graph] = make_sample_wim_graph_1();
  auto vertex_weights = [&]() {
    auto view = views::iota(vertex_id_t{10}, vertex_id_t{10} + graph::num_vertices(graph));
    return std::vector<vertex_weight_t>(view.begin(), view.end());
  }();

  auto coarsening_params = CoarseningParams{
      .neighbor_match_rule = NeighborMatchRule::HEM_P_MAX,
      .edge_weight_rule = EdgeWeightRule::SEPARATE_PRECISE,
      .seed_edge_weight_rule = SeedEdgeWeightRule::BEST_SEED_INDEX,
      .in_out_heuristic_rule = InOutHeuristicRule::P,
      .vertex_weight_rule = VertexWeightRule::AVERAGE_BY_PATHS,
  };
  auto expanding_params = ExpandingParams{
      .vertex_expanding_rule = VertexExpandingRule::ITERATIVE, .n_iterations = 5, .simulation_try_count = 5};

  auto bidir_graph = merge_wim_edge_to_undirected(graph, coarsening_params);
  ELOG_INFO << [&] {
    constexpr auto msg_pattern_header = "Merged bidirectional graph: |V|, |E| = {}";
    auto res = fmt::format(msg_pattern_header, graph_n_m(bidir_graph));
    for (auto [u, v, p] : graph::make_edge_range<0>(bidir_graph)) {
      res += fmt::format("\n\tu = {}, v = {}, p = {}", u, v, p);
    }
    return res;
  }();

  auto n_groups = 4;
  auto group_id = std::vector<vertex_id_t>{0, 1, 2, 0, 1, 2, 0, 1, 2, 3};
  // auto [n_groups, group_id] = mongoose_match(bidir_graph, coarsening_params);

  auto detailed_res = coarsen_wim_graph_by_match_d( //
      graph, inv_graph, vertex_weights, n_groups, group_id, coarsening_params);
  ELOGFMT(INFO, "Detailed coarsening result: {:4}", *detailed_res);

  auto [graph_left, inv_graph_left] = make_sample_wim_graph_1_left();
  ELOG_INFO << "Testing left part with simulation:";
  monte_carlo_test(graph_left);

  auto [graph_right, inv_graph_right] = make_sample_wim_graph_1_right();
  ELOG_INFO << "Testing right part with simulation:";
  monte_carlo_test(graph_right);

  return 0;
}

/*
Parameters: {
    "boosted_edge_weight_rule": "BEST_BOOSTED_INDEX",
    "edge_weight_rule": "SEPARATE_SIMPLE",
    "in_out_heuristic_rule": "P",
    "max_distance_from_seed": 6,
    "neighbor_match_rule": "HEM_P_MAX",
    "seed_edge_weight_rule": "BEST_SEED_INDEX",
    "vertex_weight_rule": "AVERAGE_BY_PATHS"
}
.coarsened.adj_list = {
    (0, 1): {"p":0.4022427797317505,"p_seed":0.7411531805992126},
    (1, 0): {"p":0.2200000286102295,"p_seed":0.4399999976158142},
    (1, 2): {"p":0.4385894536972046,"p_seed":0.7959796190261841},
    (2, 0): {"p":0.09320050477981567,"p_seed":0.10119998455047607},
    (2, 1): {"p":0.0,"p_seed":0.46000003814697266},
    (3, 0): {"p":0.0,"p_seed":0.9520000219345093}
},
.coarsened.inv_adj_list = {
    (0, 1): {"p":0.2200000286102295,"p_seed":0.4399999976158142},
    (0, 2): {"p":0.09320050477981567,"p_seed":0.10119998455047607},
    (0, 3): {"p":0.0,"p_seed":0.9520000219345093},
    (1, 0): {"p":0.4022427797317505,"p_seed":0.7411531805992126},
    (1, 2): {"p":0.0,"p_seed":0.46000003814697266},
    (2, 1): {"p":0.4385894536972046,"p_seed":0.7959796190261841}
},
.coarsened.vertex_weights = {
    16.5917, 19.1358, 18.8621, 19.0000
},

Parameters: {
    "boosted_edge_weight_rule": "BEST_BOOSTED_INDEX",
    "edge_weight_rule": "MERGED_PRECISE",
    "in_out_heuristic_rule": "P",
    "max_distance_from_seed": 6,
    "neighbor_match_rule": "HEM_P_MAX",
    "seed_edge_weight_rule": "BEST_SEED_INDEX",
    "vertex_weight_rule": "AVERAGE_BY_PATHS"
}
.coarsened.adj_list = {
    (0, 1): {"p":0.4730566442012787,"p_seed":0.7411531805992126},
    (1, 0): {"p":0.2200000286102295,"p_seed":0.4399999976158142},
    (1, 2): {"p":0.5040860772132874,"p_seed":0.7959796190261841},
    (2, 0): {"p":0.1167500913143158,"p_seed":0.1167500913143158},
    (2, 1): {"p":0.0,"p_seed":0.46000003814697266},
    (3, 0): {"p":0.0,"p_seed":0.9520000219345093}
},
.coarsened.inv_adj_list = {
    (0, 1): {"p":0.2200000286102295,"p_seed":0.4399999976158142},
    (0, 2): {"p":0.1167500913143158,"p_seed":0.1167500913143158},
    (0, 3): {"p":0.0,"p_seed":0.9520000219345093},
    (1, 0): {"p":0.4730566442012787,"p_seed":0.7411531805992126},
    (1, 2): {"p":0.0,"p_seed":0.46000003814697266},
    (2, 1): {"p":0.5040860772132874,"p_seed":0.7959796190261841}
},
.coarsened.vertex_weights = {
    16.5917, 19.1358, 18.8621, 19.0000
},

Parameters: {
    "boosted_edge_weight_rule": "BEST_BOOSTED_INDEX",
    "edge_weight_rule": "SEPARATE_PRECISE",
    "in_out_heuristic_rule": "P",
    "max_distance_from_seed": 6,
    "neighbor_match_rule": "HEM_P_MAX",
    "seed_edge_weight_rule": "BEST_SEED_INDEX",
    "vertex_weight_rule": "AVERAGE_BY_PATHS"
}
.coarsened.adj_list = {
    (0, 1): {"p":0.39859136939048767,"p_seed":0.7411531805992126},
    (1, 0): {"p":0.2200000286102295,"p_seed":0.4399999976158142},
    (1, 2): {"p":0.4350780248641968,"p_seed":0.7959796190261841},
    (2, 0): {"p":0.09320050477981567,"p_seed":0.10119998455047607},
    (2, 1): {"p":0.0,"p_seed":0.46000003814697266},
    (3, 0): {"p":0.0,"p_seed":0.9520000219345093}
},
.coarsened.inv_adj_list = {
    (0, 1): {"p":0.2200000286102295,"p_seed":0.4399999976158142},
    (0, 2): {"p":0.09320050477981567,"p_seed":0.10119998455047607},
    (0, 3): {"p":0.0,"p_seed":0.9520000219345093},
    (1, 0): {"p":0.39859136939048767,"p_seed":0.7411531805992126},
    (1, 2): {"p":0.0,"p_seed":0.46000003814697266},
    (2, 1): {"p":0.4350780248641968,"p_seed":0.7959796190261841}
},
.coarsened.vertex_weights = {
    16.5917, 19.1358, 18.8621, 19.0000
},

Testing left part with simulation:
    seeds = [0], 1, dest = [1, 3, 5], 101010
    Simulated result with seed 0 = 0.421093
    seeds = [2], 100, dest = [1, 3, 5], 101010
    Simulated result with seed 2 = 0.421757
    seeds = [4], 10000, dest = [1, 3, 5], 101010
    Simulated result with seed 4 = 0.339945
    seeds = [6], 1000000, dest = [1, 3, 5], 101010
    Simulated result with seed 6 = 0.336341
    seeds = [0, 2], 101, dest = [1, 3, 5], 101010
    Simulated result with seed [0, 2] = 0.595138
    seeds = [0, 4], 10001, dest = [1, 3, 5], 101010
    Simulated result with seed [0, 4] = 0.598495
    seeds = [2, 4], 10100, dest = [1, 3, 5], 101010

Testing right part with simulation:
    seeds = [0], 1, dest = [1, 3, 5], 101010
    Simulated result with seed 0 = 0.435929
    seeds = [2], 100, dest = [1, 3, 5], 101010
    Simulated result with seed 2 = 0.470347
    seeds = [4], 10000, dest = [1, 3, 5], 101010
    Simulated result with seed 4 = 0.400862
    seeds = [6], 1000000, dest = [1, 3, 5], 101010
    Simulated result with seed 6 = 0.354797
    seeds = [0, 2], 101, dest = [1, 3, 5], 101010
    Simulated result with seed [0, 2] = 0.610634
    seeds = [0, 4], 10001, dest = [1, 3, 5], 101010
    Simulated result with seed [0, 4] = 0.631561
    seeds = [2, 4], 10100, dest = [1, 3, 5], 101010
    Simulated result with seed [2, 4] = 0.591601
*/
