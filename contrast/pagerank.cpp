#include "contrast_algorithms.h"
#include "utils/easylog.h"

namespace {
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

        takes_max ? ranges::stable_sort(indices, ranges::greater{}, LAMBDA_1(PR[_1]))
                  : ranges::stable_sort(indices, ranges::less{}, LAMBDA_1(PR[_1]));
        return std::vector<vertex_id_t>{indices.begin(), indices.begin() + k};
      });
}
} // namespace

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
