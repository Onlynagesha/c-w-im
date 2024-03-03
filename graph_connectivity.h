#pragma once

#include "utils/dynamic_bitset.h"
#include "utils/graph.h"
#include "utils/utils.h"
#include <queue>
#include <stack>

namespace details {
// Counts with Tarjan's SCC algorithm
struct SCCTarjanContext {
  static constexpr auto NULL_INDEX = std::numeric_limits<vertex_id_t>::max();

  vertex_id_t n;
  vertex_id_t result;
  vertex_id_t next_timestamp;
  DynamicBitset in_stack;
  std::stack<vertex_id_t> stack;
  std::vector<vertex_id_t> dfs_timestamp;
  std::vector<vertex_id_t> low;

  explicit SCCTarjanContext(vertex_id_t n)
      : n(n), result(0), next_timestamp(0), in_stack(n), dfs_timestamp(n, NULL_INDEX), low(n, NULL_INDEX) {}
};

template <class... Attributes>
inline auto n_strongly_connected_components_dfs(SCCTarjanContext& ctx, const AdjacencyList<Attributes...>& graph,
                                                vertex_id_t v) -> void {
  ctx.dfs_timestamp[v] = ctx.low[v] = ctx.next_timestamp++;
  ctx.in_stack.set(v);
  ctx.stack.push(v);
  for (auto u : graph[v] | views::keys) {
    if (ctx.dfs_timestamp[u] == SCCTarjanContext::NULL_INDEX) {
      n_strongly_connected_components_dfs(ctx, graph, u);
      ctx.low[v] = std::min(ctx.low[v], ctx.low[u]);
    } else if (ctx.in_stack.test(u)) {
      ctx.low[v] = std::min(ctx.low[v], ctx.dfs_timestamp[u]); // Back edge detected
    }
  }
  // Counts for each DFS root
  if (ctx.low[v] == ctx.dfs_timestamp[v]) {
    ctx.result += 1;
    while (true) {
      BOOST_ASSERT(!ctx.stack.empty());
      auto top = ctx.stack.top();
      ctx.stack.pop();
      ctx.in_stack.reset(top);
      if (top == v) {
        break;
      }
    }
  }
}
} // namespace details

template <class... Attributes>
inline auto n_strongly_connected_components(const AdjacencyList<Attributes...>& graph) {
  auto ctx = details::SCCTarjanContext(graph::num_vertices(graph));
  for (auto v : vertices(graph)) {
    if (ctx.dfs_timestamp[v] == details::SCCTarjanContext::NULL_INDEX) {
      details::n_strongly_connected_components_dfs(ctx, graph, v);
    }
  }
  return ctx.result;
}

template <class... Attributes>
inline auto n_weakly_connected_components(const AdjacencyList<Attributes...>& graph,
                                          const InvAdjacencyList<Attributes...>& inv_graph) {
  auto n = graph::num_vertices(graph);
  auto visited = DynamicBitset(n);
  auto queue = std::queue<vertex_id_t>{};

  auto push_to_queue = [&](vertex_id_t v) {
    queue.push(v);
    visited.set(v);
  };
  auto pop_from_queue = [&]() {
    auto res = queue.front();
    queue.pop();
    return res;
  };

  auto res = vertex_id_t{0};
  for (auto v : vertices(graph)) {
    if (visited.test(v)) {
      continue;
    }
    // Starts BFS from v
    push_to_queue(v);
    while (!queue.empty()) {
      auto cur = pop_from_queue();
      for (auto out_neighbor : graph[cur] | views::keys) {
        if (!visited.test(out_neighbor)) {
          push_to_queue(out_neighbor);
        }
      }
      for (auto in_neighbor : inv_graph[cur] | views::keys) {
        if (!visited.test(in_neighbor)) {
          push_to_queue(in_neighbor);
        }
      }
    }
    res += 1; // New component {v, ...} counted
  }
  return res;
}