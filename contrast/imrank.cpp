#include "contrast_algorithms.h"
#include "utils/easylog.h"
#include <fmt/ranges.h>

namespace {
template <is_edge_property E>
auto imrank_generic(const InvAdjacencyList<E>& inv_graph, const IMRankParams& params) -> std::vector<vertex_id_t> {
  auto n = graph::num_vertices(inv_graph);
  auto k = (params.k <= 0) ? n : std::min(params.k, n);
  // In descending order of rank
  auto sorted_vertices = max_out_degree(inv_graph);
  auto sorted_vertices_temp = std::vector<vertex_id_t>(n);
  MYLOG_FMT_TRACE("Initial sorted_vertices = {}", sorted_vertices);
  // The lower rank[v] value, the higher rank v is
  auto rank = std::vector<vertex_id_t>(n);
  // Marginal influence spread
  auto Mr = std::vector<edge_probability_t>();

  auto topk_unchanged_count = uint64_t{0};
  for (auto iteration_index = uint64_t{0}; iteration_index < params.n_iterations; ++iteration_index) {
    Mr.assign(n, 1.0_ep);
    // Resets rank of each vertex
    for (auto [rank_v, v] : sorted_vertices | views::enumerate) {
      rank[v] = rank_v;
    }
    MYLOG_FMT_TRACE("rank => {}", rank);
    // From the lowest rank to the highest rank
    for (auto [rank_vi, vi] : sorted_vertices | views::enumerate | views::reverse) {
      // Traverses each in-edge (vj, vi) such that rank of vj is higher than rank of vi
      for (auto [vj, w] : inv_graph[vi]) {
        if (rank[vj] >= rank_vi) {
          continue; // Skips the in-neighbors with lower rank
        }
        auto delta = w.p * Mr[vi];
        Mr[vj] += delta;
        Mr[vi] -= delta;
      }
    }
    MYLOG_FMT_TRACE("At iteration #{}: Mr = {::.4f}", iteration_index, Mr);

    ranges::copy(sorted_vertices, sorted_vertices_temp.begin());
    ranges::stable_sort(sorted_vertices_temp, ranges::greater{}, LAMBDA_1(Mr[_1]));
    MYLOG_FMT_TRACE("At iteration #{}: sorted_vertices => {}", iteration_index, sorted_vertices_temp);

    if (ranges::equal(sorted_vertices | views::take(k), sorted_vertices_temp | views::take(k))) {
      topk_unchanged_count += 1;
      MYLOG_FMT_DEBUG("At iteration #{}: Top-{} keeps unchanged in {} iterations.", //
                      iteration_index, k, topk_unchanged_count);
    } else {
      topk_unchanged_count = 0;
      MYLOG_FMT_TRACE("At iteration #{}: Top-{} changes to {}", //
                      iteration_index, k, sorted_vertices_temp | views::take(k));
    }
    // sorted_vertices <- sorted_vertices_temp
    sorted_vertices.swap(sorted_vertices_temp);
    // Early stop of top-k keeps unchanged for sufficient iterations
    if (topk_unchanged_count >= params.n_iterations_before_topk_fixed) {
      MYLOG_FMT_DEBUG("Early stop after {} iterations.", iteration_index + 1);
      break;
    }
  }
  // Takes the top-k
  return {sorted_vertices.begin(), sorted_vertices.begin() + k};
}
} // namespace

auto wim_imrank(const InvAdjacencyList<WIMEdge>& inv_graph, const IMRankParams& params) -> std::vector<vertex_id_t> {
  return imrank_generic(inv_graph, params);
}

auto wbim_imrank(const InvAdjacencyList<WBIMEdge>& inv_graph, const IMRankParams& params) -> std::vector<vertex_id_t> {
  return imrank_generic(inv_graph, params);
}
