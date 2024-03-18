#include "contrast_algorithms.h"
#include "utils/easylog.h"
#include <fmt/ranges.h>

namespace {
template <is_edge_property E, int IsInv>
auto max_out_degree_generic(const graph::adjacency<IsInv, E>& graph, vertex_id_t top_k, bool uses_edge_weight)
    -> std::vector<vertex_id_t> {
  auto n = graph::num_vertices(graph);
  if (top_k <= 0) {
    top_k = n;
  } else {
    top_k = std::min(top_k, n);
  }
  auto out_degrees = [&]() {
    auto res = std::vector<edge_probability_t>(n, 0.0_ep);
    if constexpr (IsInv == 0) {
      // Traverses out-edges for each vertex v
      for (auto source : vertices(graph)) {
        for (auto [target, w] : graph[source]) {
          res[source] += (uses_edge_weight ? w.p : 1.0_ep);
        }
      }
    } else {
      // Otherwise, traverses each in-edge in the transpose graph
      for (auto target : vertices(graph)) {
        for (auto [source, w] : graph[target]) {
          res[source] += (uses_edge_weight ? w.p : 1.0_ep);
        }
      }
    }
    return res;
  }();
  auto indices = std::vector<vertex_id_t>(n);
  std::iota(indices.begin(), indices.end(), 0);
  // Sorts indices by (weighted) out-degree in descending order
  ranges::sort(indices, ranges::greater{}, LAMBDA_1(out_degrees[_1]));
  MYLOG_TRACE([&]() {
    auto str = fmt::format("Result of {} sorting:", uses_edge_weight ? "max-strength" : "max-degree");
    auto arr_name = uses_edge_weight ? "max_strength" : "max_degree";
    for (auto i : indices) {
      str += fmt::format("\n\t{}[{}] = {:.4f}", arr_name, i, out_degrees[i]);
    }
    return str;
  }());
  // Takes the top-k
  return {indices.begin(), indices.begin() + top_k};
}

template <is_edge_property E>
auto pagerank_generic(const AdjacencyList<E>& graph, const InvAdjacencyList<E>& inv_graph,
                      std::span<const edge_probability_t> vertex_weights, const PagerankParams& params)
    -> rfl::Result<std::vector<edge_probability_t>> {
  auto n = graph::num_vertices(graph);
  if (n <= 0) {
    return std::vector<edge_probability_t>{}; // Returns empty list for empty graph
  }
  // Checking n and vertex_weights
  if (params.uses_vertex_weight) {
    if (vertex_weights.size() != n) {
      constexpr auto msg_pattern = "Size mismatch between |V| (which is {}) "
                                   "and the length of vertex weights list (which is {})";
      return rfl::Error{fmt::format(msg_pattern, n, vertex_weights.size())};
    }
    if (auto min_vertex_weight = ranges::min(vertex_weights);
        std::isnan(min_vertex_weight) || std::isinf(min_vertex_weight) || min_vertex_weight <= 0.0) {
      constexpr auto msg_pattern = "Non-positive vertex weight (which is {}) is not supported.";
      return rfl::Error{fmt::format(msg_pattern, min_vertex_weight)};
    }
  }

  constexpr auto to_square = views::transform([](auto x) { return x * x; });
  constexpr auto vector_l2_norm = [to_square](std::span<const edge_probability_t> values) {
    auto sqr_sum = accumulate_sum(values | to_square);
    BOOST_ASSERT_MSG(!std::isnan(sqr_sum) && !std::isinf(sqr_sum) && sqr_sum > 0.0, "Expects ||x|| > 0.");
    return std::sqrt(sqr_sum);
  };

  // Initially PR(v) = c(v) for each v.
  auto [res, sum_vertex_weights] = [&]() {
    if (params.uses_vertex_weight) {
      // c(v) = vertex_weights[v] for each vertex
      auto res = std::vector(vertex_weights.begin(), vertex_weights.end());
      return std::tuple{std::move(res), accumulate_sum(vertex_weights)};
    }
    // Otherwise, c(v) 1.0 for each vertex
    return std::tuple{std::vector<edge_probability_t>(n, 1.0), static_cast<edge_probability_t>(n)};
  }();
  auto res_l2_norm = vector_l2_norm(res);
  MYLOG_FMT_TRACE("PR(0) = {::.4f}, with L2-norm = {:.4f}", res, res_l2_norm);
  // L(v) = Sum of (v -> u).p (if weighted) for each out-neighbor u if transpose == false
  // L(v) = Sum of (x -> v).p (if weighted) for each in-neighbor x if transpose == true
  auto total_d = [&]() {
    auto get_fn = [&params, n](const auto& which_graph) {
      auto res = std::vector<edge_probability_t>(n, 0.0_ep);
      for (auto u : vertices(which_graph)) {
        for (auto [v, w] : which_graph[u]) {
          res[u] += (params.uses_edge_weight ? w.p : 1.0_ep);
        }
      }
      return res;
    };
    return params.transpose ? get_fn(inv_graph) : get_fn(graph);
  }();

  auto temp = std::vector<edge_probability_t>(n);
  // PR'(v) <- Sum of PR(u) * (u -> v).p / L(u) for each in-neighbor u if edge weight is used
  // PR'(v) <- Sum of PR(u) / L(u) for each in-neighbor u otherwise
  for (auto iteration_index = uint64_t{0}; iteration_index < params.n_iterations; iteration_index++) {
    for (auto v : range(n)) {
      auto sum = 0.0;
      auto loop_body = [&](auto&& edge) {
        auto [u, w] = edge;
        if (params.uses_edge_weight) {
          if (w.p > 0) {
            BOOST_ASSERT(0 < w.p && w.p <= total_d[u]);
            sum += res[u] * w.p / total_d[u];
          }
        } else {
          BOOST_ASSERT(0 < total_d[u]);
          sum += res[u] / total_d[u];
        }
      };
      params.transpose ? ranges::for_each(graph[v], loop_body) : ranges::for_each(inv_graph[v], loop_body);
      // PR(v) <- d * PR(v) + (1 - d) * c(v)
      temp[v] = *params.damping_factor * sum +
                (1.0_ep - *params.damping_factor) * (params.uses_vertex_weight ? vertex_weights[v] : 1.0_ep);
    }
    // Compares changed by L2-norm
    auto temp_l2_norm = vector_l2_norm(temp);
    auto early_stops = [&]() {
      auto diff_sqr_l2_norm = //  || x0 - x1 ||^2 where x0 = res / ||res||, x1 = temp / ||temp||
          accumulate_sum(range(n) | TRANSFORM_VIEW(res[_1] / res_l2_norm - temp[_1] / temp_l2_norm) | to_square);
      MYLOG_FMT_TRACE("L2-norm of PR({}) = {:.4f}; Squared L2-norm of difference = {:.4e}", //
                      iteration_index + 1, temp_l2_norm, diff_sqr_l2_norm);
      return diff_sqr_l2_norm < (*params.epsilon) * (*params.epsilon);
    }();
    // res <- temp
    res.swap(temp);
    res_l2_norm = temp_l2_norm;
    MYLOG_FMT_TRACE("PR({}) = {::.4f}", iteration_index + 1, res);
    if (early_stops) {
      MYLOG_FMT_DEBUG("Pagerank early-stops after {} iterations due to convergence.", iteration_index + 1);
      break;
    }
  }
  return res;
}

template <is_edge_property E>
auto best_pagerank_generic(const AdjacencyList<E>& graph, const InvAdjacencyList<E>& inv_graph,
                           std::span<const edge_probability_t> vertex_weights, const PagerankParams& params,
                           bool takes_max) -> rfl::Result<std::vector<vertex_id_t>> {
  auto n = graph::num_vertices(graph);
  auto k = (params.k <= 0) ? n : std::min(params.k, n);
  return pagerank_generic(graph, inv_graph, vertex_weights, params)
      .transform([n, k, takes_max](std::vector<edge_probability_t> PR) {
        auto indices = std::vector<vertex_id_t>(n);
        std::iota(indices.begin(), indices.end(), 0);

        takes_max ? ranges::sort(indices, ranges::greater{}, LAMBDA_1(PR[_1]))
                  : ranges::sort(indices, ranges::less{}, LAMBDA_1(PR[_1]));
        return std::vector<vertex_id_t>{indices.begin(), indices.begin() + k};
      });
}

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
    ranges::sort(sorted_vertices_temp, ranges::greater{}, LAMBDA_1(Mr[_1]));
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

auto wim_max_out_degree(const AdjacencyList<WIMEdge>& graph, vertex_id_t top_k) -> std::vector<vertex_id_t> {
  return max_out_degree_generic(graph, top_k, false);
}

auto wbim_max_out_degree(const AdjacencyList<WBIMEdge>& graph, vertex_id_t top_k) -> std::vector<vertex_id_t> {
  return max_out_degree_generic(graph, top_k, false);
}

auto wim_max_out_degree_i(const InvAdjacencyList<WIMEdge>& inv_graph, vertex_id_t top_k) -> std::vector<vertex_id_t> {
  return max_out_degree_generic(inv_graph, top_k, false);
}

auto wbim_max_out_degree_i(const InvAdjacencyList<WBIMEdge>& inv_graph, vertex_id_t top_k) -> std::vector<vertex_id_t> {
  return max_out_degree_generic(inv_graph, top_k, false);
}

auto wim_max_out_strength(const AdjacencyList<WIMEdge>& graph, vertex_id_t top_k) -> std::vector<vertex_id_t> {
  return max_out_degree_generic(graph, top_k, true);
}

auto wbim_max_out_strength(const AdjacencyList<WBIMEdge>& graph, vertex_id_t top_k) -> std::vector<vertex_id_t> {
  return max_out_degree_generic(graph, top_k, true);
}

auto wim_max_out_strength_i(const InvAdjacencyList<WIMEdge>& inv_graph, vertex_id_t top_k) -> std::vector<vertex_id_t> {
  return max_out_degree_generic(inv_graph, top_k, true);
}

auto wbim_max_out_strength_i(const InvAdjacencyList<WBIMEdge>& inv_graph, vertex_id_t top_k)
    -> std::vector<vertex_id_t> {
  return max_out_degree_generic(inv_graph, top_k, true);
}

auto wim_pagerank(const AdjacencyList<WIMEdge>& graph, const InvAdjacencyList<WIMEdge>& inv_graph,
                  std::span<const edge_probability_t> vertex_weights, const PagerankParams& params)
    -> rfl::Result<std::vector<edge_probability_t>> {
  return pagerank_generic(graph, inv_graph, vertex_weights, params);
}

auto wbim_pagerank(const AdjacencyList<WBIMEdge>& graph, const InvAdjacencyList<WBIMEdge>& inv_graph,
                   std::span<const edge_probability_t> vertex_weights, const PagerankParams& params)
    -> rfl::Result<std::vector<edge_probability_t>> {
  return pagerank_generic(graph, inv_graph, vertex_weights, params);
}

auto wim_min_pagerank(const AdjacencyList<WIMEdge>& graph, const InvAdjacencyList<WIMEdge>& inv_graph,
                      std::span<const edge_probability_t> vertex_weights, const PagerankParams& params)
    -> rfl::Result<std::vector<vertex_id_t>> {
  return best_pagerank_generic(graph, inv_graph, vertex_weights, params, false);
}

auto wim_max_pagerank(const AdjacencyList<WIMEdge>& graph, const InvAdjacencyList<WIMEdge>& inv_graph,
                      std::span<const edge_probability_t> vertex_weights, const PagerankParams& params)
    -> rfl::Result<std::vector<vertex_id_t>> {
  return best_pagerank_generic(graph, inv_graph, vertex_weights, params, true);
}

auto wbim_min_pagerank(const AdjacencyList<WBIMEdge>& graph, const InvAdjacencyList<WBIMEdge>& inv_graph,
                       std::span<const edge_probability_t> vertex_weights, const PagerankParams& params)
    -> rfl::Result<std::vector<vertex_id_t>> {
  return best_pagerank_generic(graph, inv_graph, vertex_weights, params, false);
}

auto wbim_max_pagerank(const AdjacencyList<WBIMEdge>& graph, const InvAdjacencyList<WBIMEdge>& inv_graph,
                       std::span<const edge_probability_t> vertex_weights, const PagerankParams& params)
    -> rfl::Result<std::vector<vertex_id_t>> {
  return best_pagerank_generic(graph, inv_graph, vertex_weights, params, true);
}

auto wim_imrank(const InvAdjacencyList<WIMEdge>& inv_graph, const IMRankParams& params) -> std::vector<vertex_id_t> {
  return imrank_generic(inv_graph, params);
}

auto wbim_imrank(const InvAdjacencyList<WBIMEdge>& inv_graph, const IMRankParams& params) -> std::vector<vertex_id_t> {
  return imrank_generic(inv_graph, params);
}
