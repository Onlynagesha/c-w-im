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
  ranges::stable_sort(indices, ranges::greater{}, LAMBDA_1(out_degrees[_1]));
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