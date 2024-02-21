#include "coarsening.h"
#include "dump.h"
#include "playground/sample_graph.h"
#include "utils/graph_connectivity.h"
#include <fmt/ranges.h>
#include <nwgraph/adaptors/edge_range.hpp>
#include <nwgraph/adaptors/neighbor_range.hpp>
#include <ylt/easylog.hpp>

int main() {
  auto irr_graph = DirectedEdgeList<double>{};
  irr_graph.open_for_push_back();
  irr_graph.push_back(5, 10, 1.25);
  irr_graph.push_back(10, 5, 1.75);
  irr_graph.push_back(10, 20, 2.25);
  irr_graph.close_for_push_back();

  auto irr_adj_list = AdjacencyList<double>{irr_graph};
  ELOGFMT(INFO, "irr_adj_list: |V| = {}, |E| = {}", irr_adj_list.num_vertices()[0], irr_adj_list.num_edges());
  for (auto [u, v, w] : nw::graph::make_edge_range<0>(irr_adj_list)) {
    ELOGFMT(INFO, "\t({}, {}): {}", u, v, w);
  }
  auto irr_n_scc = n_strongly_connected_components(irr_adj_list);
  auto irr_n_wcc = n_weakly_connected_components(irr_adj_list, InvAdjacencyList<double>{irr_graph});
  ELOGFMT(INFO, "# of SCC = {}; # of WCC = {}", irr_n_scc, irr_n_wcc);

  auto [graph, inv_graph] = make_sample_wim_graph();
  auto coarsening_params = CoarseningParams{.neighbor_match_rule = NeighborMatchRule::HEM_P_MAX,
                                            .group_path_probability_rule = GroupPathProbabilityRule::P_SEPARATE,
                                            .group_in_out_rule = GroupInOutRule::W,
                                            .vertex_weight_rule = VertexWeightRule::AVERAGE_BY_PATHS,
                                            .seed_merge_rule = SeedMergeRule::S_SINGLE};
  auto expanding_params = ExpandingParams{
      .seed_expansion_rule = SeedExpansionRule::S_ITERATIVE, .n_iterations = 10, .simulation_try_count = 10};

  auto bidir_graph = merge_wim_edge_to_undirected(graph, coarsening_params);
  ELOG_INFO << [&] {
    constexpr auto msg_pattern_header = "Merged bidirectional graph (|V| = {}, |E| = {}):";
    auto res = fmt::format(msg_pattern_header, bidir_graph.num_vertices()[0], bidir_graph.num_edges());
    for (auto [u, v, p] : graph::make_edge_range<0>(bidir_graph)) {
      res += fmt::format("\n\tu = {}, v = {}, p = {}", u, v, p);
    }
    return res;
  }();

  // Without seeds
  {
    auto group_result = mongoose_match(bidir_graph, coarsening_params);
    ELOG_INFO << [&] {
      auto res = fmt::format("Group result #1 (with n_groups = {}):", group_result.n_groups);
      for (auto [v, g] : views::enumerate(group_result.group_id)) {
        res += fmt::format("\n\t{}: {}", v, g);
      }
      return res;
    }();
    auto vertex_weights = [&]() {
      auto view = views::iota(vertex_id_t{10}, vertex_id_t{10} + graph::num_vertices(graph));
      return std::vector<vertex_weight_t>(view.begin(), view.end());
    }();
    auto coarsening_res = do_coarsen_wim_graph_w(graph, inv_graph, vertex_weights, group_result.n_groups,
                                                 group_result.group_id, coarsening_params);
    ELOGFMT(INFO, "Coarsening result: {:4}", coarsening_res);
  }
  ELOG_INFO << "==============================================================";
  // With seeds
  {
    auto seeds = {1_vid, 7_vid};
    auto group_result = mongoose_match_with_seeds(bidir_graph, seeds, coarsening_params);
    ELOG_INFO << [&] {
      auto res = fmt::format("Group result #2 (with n_groups = {}, seeds = {}):", group_result.n_groups, seeds);
      for (auto [v, g] : views::enumerate(group_result.group_id)) {
        res += fmt::format("\n\t{}: {}", v, g);
      }
      return res;
    }();
    auto vertex_weights = [&]() {
      auto view = views::iota(vertex_id_t{10}, vertex_id_t{10} + graph::num_vertices(graph));
      return std::vector<vertex_weight_t>(view.begin(), view.end());
    }();
    auto coarsening_res = do_coarsen_wim_graph_w(graph, inv_graph, vertex_weights, group_result.n_groups,
                                                 group_result.group_id, coarsening_params);
    ELOGFMT(INFO, "Coarsening result: {:4}", coarsening_res);

    auto coarsened_seeds = std::vector<vertex_id_t>{0, 2, 3};

    expand_wim_seed_vertices_w(graph, vertex_weights, coarsening_res, coarsened_seeds, expanding_params)
        .and_then([&](ExpandSeedResult expand_res) -> ResultVoid {
          ELOGFMT(INFO, "Expanded seeds: {}", expand_res.expanded_seeds);
          return RESULT_VOID_SUCCESS;
        });
  }

  return 0;
}
