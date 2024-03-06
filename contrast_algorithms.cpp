#include "contrast_algorithms.h"
#include "utils/easylog.h"
#include <fmt/ranges.h>

auto wim_pagerank(const AdjacencyList<WIMEdge>& graph, const InvAdjacencyList<WIMEdge>& inv_graph,
                  std::span<const edge_probability_t> vertex_weights, const PagerankParams& params)
    -> rfl::Result<std::vector<edge_probability_t>> {
  auto n = graph::num_vertices(graph);
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
  MYLOG_FMT_TRACE("PR(0) = {::.4f}", res);
  // L(v) = (weighted) out-degree of v
  auto total_d_out = [&]() {
    if (params.uses_edge_weight) {
      // L(v) = Sum of (v -> u).p for each out-neighbor u
      auto view = graph | TRANSFORM_VIEW(accumulate_sum(_1 | TRANSFORM_VIEW(get<1>(_1).p)));
      return std::vector(view.begin(), view.end());
    }
    // Otherwise, L(v) = out-degree(v)
    auto unweighted_view = graph | TRANSFORM_VIEW(static_cast<edge_probability_t>(ranges::size(_1)));
    return std::vector(unweighted_view.begin(), unweighted_view.end());
  }();

  auto temp = std::vector<edge_probability_t>(n);
  for (auto iteration_index = uint64_t{0}; iteration_index < params.n_iterations; iteration_index++) {
    for (auto v : range(n)) {
      // PR'(v) <- Sum of PR(u) * (u -> v).p / L(u) for each in-neighbor u if edge weight is used
      // PR'(v) <- Sum of PR(u) / L(u) for each in-neighbor u otherwise
      auto sum = 0.0;
      for (auto [u, w] : inv_graph[v]) {
        sum += res[u] * (params.uses_edge_weight ? w.p : 1.0_ep) / total_d_out[u];
      }
      // PR(v) <- d * PR(v) + (1 - d) * c(v)
      temp[v] = *params.damping_factor * sum +
                (1.0_ep - *params.damping_factor) * (params.uses_vertex_weight ? vertex_weights[v] : 1.0_ep);
    }
    res.swap(temp); // res <- temp
    MYLOG_FMT_TRACE("PR({}) = {::.4f}", iteration_index + 1, res);
  }
  return res;
}
