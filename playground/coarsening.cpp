#include "coarsening.h"
#include "dump.h"
#include "graph_connectivity.h"
#include "playground/sample_graph.h"
#include "utils/easylog.h"
#include "utils/histogram.h"
#include "wim.h"
#include <fmt/ranges.h>
#include <nwgraph/adaptors/edge_range.hpp>
#include <nwgraph/adaptors/neighbor_range.hpp>

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
    for (auto [u, v, p] : graph::make_edge_range<0>(remove_const(graph))) {
      res += fmt::format("\n\tu = {}, v = {}, p = {}", u, v, p);
    }
    return res;
  }();
  for (auto seed : std::vector<vertex_id_t>{0, 2, 4, 6}) {
    auto result = monte_carlo_simulate(graph, {7, {seed}}, {7, {1, 3, 5}}, 2'000'000uLL);
    ELOGFMT(INFO, "Simulated result with seed {} = {:.6f}", seed, result);
  }
  for (auto seeds : std::vector<VertexSet>{{7, {0, 2}}, {7, {0, 4}}, {7, {2, 4}}}) {
    auto result = monte_carlo_simulate(graph, seeds, {7, {1, 3, 5}}, 2'000'000uLL);
    ELOGFMT(INFO, "Simulated result with seed {} = {:.6f}", seeds.vertex_list, result);
  }
}

int main() {
  easylog::set_min_severity(easylog::Severity::TRACE);

  auto rand_values = range(1'000'000) | TRANSFORM_VIEW(rand_float());
  ELOGFMT(INFO, "Histogram of rand_values:\n{}", make_histogram(rand_values, 100, 20));

  auto [graph, inv_graph] = make_sample_wim_graph_1();
  auto coarsening_params = CoarseningParams{.neighbor_match_rule = NeighborMatchRule::LEM_P_PRODUCT,
                                            .edge_weight_rule = EdgeWeightRule::SEPARATE_SIMPLE,
                                            .edge_seed_weight_rule = EdgeSeedWeightRule::BEST_SEED_INDEX,
                                            .in_out_heuristic_rule = InOutHeuristicRule::COUNT,
                                            .vertex_weight_rule = VertexWeightRule::AVERAGE,
                                            .seed_merging_rule = SeedMergingRule::UNUSED};
  auto expanding_params = ExpandingParams{
      .seed_expanding_rule = SeedExpandingRule::ITERATIVE, .n_iterations = 5, .simulation_try_count = 5};

  auto bidir_graph = merge_wim_edge_to_undirected(graph, coarsening_params);
  ELOG_INFO << [&] {
    constexpr auto msg_pattern_header = "Merged bidirectional graph: |V|, |E| = {}";
    auto res = fmt::format(msg_pattern_header, graph_n_m(bidir_graph));
    for (auto [u, v, p] : graph::make_edge_range<0>(bidir_graph)) {
      res += fmt::format("\n\tu = {}, v = {}, p = {}", u, v, p);
    }
    return res;
  }();

  auto vertex_weights = [&]() {
    auto view = views::iota(vertex_id_t{10}, vertex_id_t{10} + graph::num_vertices(graph));
    return std::vector<vertex_weight_t>(view.begin(), view.end());
  }();
  auto n_groups = 4;
  auto group_id = std::vector<vertex_id_t>{0, 1, 2, 0, 1, 2, 0, 1, 2, 3};
  // auto [n_groups, group_id] = mongoose_match(bidir_graph, coarsening_params);

  auto brief_res = coarsen_wim_graph_with_match_result_w( //
      graph, inv_graph, vertex_weights, n_groups, group_id, coarsening_params);
  ELOGFMT(INFO, "Brief coarsening result: {:4}", brief_res);

  auto detailed_res = coarsen_wim_graph_with_match_result_d_w( //
      graph, inv_graph, vertex_weights, n_groups, group_id, coarsening_params);
  ELOGFMT(INFO, "Detailed coarsening result: {:4}", detailed_res);

  auto [graph_left, inv_graph_left] = make_sample_wim_graph_1_left();
  ELOG_INFO << "Testing left part with simulation:";
  monte_carlo_test(graph_left);

  auto [graph_right, inv_graph_right] = make_sample_wim_graph_1_right();
  ELOG_INFO << "Testing right part with simulation:";
  monte_carlo_test(graph_right);

  return 0;
}
