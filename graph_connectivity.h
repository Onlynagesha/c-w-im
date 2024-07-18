#pragma once

#include "utils/dynamic_bitset.h"
#include "utils/graph.h"
#include "utils/utils.h"
#include <queue>
#include <stack>

namespace details {
// Counts with Tarjan's SCC algorithm
struct SCCTarjanContext {
  static_assert(std::unsigned_integral<vertex_id_t>, "Assumes vertex_id_t to be unsigned.");
  static constexpr auto NULL_INDEX = std::numeric_limits<vertex_id_t>::max();

  struct DFSState {
    vertex_id_t v;
    vertex_id_t index;

    DFSState(vertex_id_t v) : v(v), index(NULL_INDEX) {}
  };

  // # of vertices, i.e. |V|
  vertex_id_t n;
  // # of SCCs found
  vertex_id_t result;
  // Accumulator of timestamp, i.e. # of vertices traversed before
  vertex_id_t next_timestamp;
  // in_trace_stack[v] = Whether v is in trace_stack
  DynamicBitset in_trace_stack;
  // The stack used to simulate the DFS process in a non-recursive manner.
  std::stack<DFSState> dfs_stack;
  // The stack used to trace SCCs.
  std::stack<vertex_id_t> trace_stack;
  // dfs_timestamp[v] = Time when v is traversed. Also known as dfn[v].
  std::vector<vertex_id_t> dfs_timestamp;
  // See documentation of Tarjan algorithm.
  std::vector<vertex_id_t> low;

  explicit SCCTarjanContext(vertex_id_t n)
      : n(n), result(0), next_timestamp(0), in_trace_stack(n), dfs_timestamp(n, NULL_INDEX), low(n, NULL_INDEX) {}
};

template <class... Attributes>
inline auto n_strongly_connected_components_dfs( //
    SCCTarjanContext& ctx, const AdjacencyList<Attributes...>& graph, vertex_id_t v) -> void {
  // (a) v is visited for the first time.
  ctx.dfs_timestamp[v] = ctx.low[v] = ctx.next_timestamp++;
  ctx.in_trace_stack.set(v);
  ctx.trace_stack.push(v);
  // (b) DFS
  for (auto u : graph[v] | views::keys) {
    if (ctx.dfs_timestamp[u] == SCCTarjanContext::NULL_INDEX) {
      // (b.1) Recursion
      n_strongly_connected_components_dfs(ctx, graph, u);
      // (b.2) Updates low[v] by low[u]
      ctx.low[v] = std::min(ctx.low[v], ctx.low[u]);
    } else if (ctx.in_trace_stack.test(u)) {
      // (b.3) Updates low[v] by dfn[u] on back edge detected
      ctx.low[v] = std::min(ctx.low[v], ctx.dfs_timestamp[u]);
    }
  }
  // (c) Counts for each DFS root
  if (ctx.low[v] == ctx.dfs_timestamp[v]) {
    ctx.result += 1;
    while (true) {
      BOOST_ASSERT(!ctx.stack.empty());
      auto top = ctx.trace_stack.top();
      ctx.trace_stack.pop();
      ctx.in_trace_stack.reset(top);
      if (top == v) {
        break;
      }
    }
  }
  // (d) Returns
}

template <class... Attributes>
inline auto n_strongly_connected_components_dfs_non_recursive( //
    SCCTarjanContext& ctx, const AdjacencyList<Attributes...>& graph, vertex_id_t v) -> void {
  ctx.dfs_stack.emplace(v);
  while (!ctx.dfs_stack.empty()) {
  n_strongly_connected_components_dfs_non_recursive_head:
    auto& [v, index] = ctx.dfs_stack.top();
    auto adj_list_v = graph[v] | views::keys;
    auto adj_list_size = ranges::size(adj_list_v);

    if (index == SCCTarjanContext::NULL_INDEX) {
      // (a) See above
      ctx.dfs_timestamp[v] = ctx.low[v] = ctx.next_timestamp++;
      ctx.in_trace_stack.set(v);
      ctx.trace_stack.push(v);
    } else {
      // (b.2) See above
      auto u = adj_list_v[index];
      ctx.low[v] = std::min(ctx.low[v], ctx.low[u]);
    }
    while ((++index) < adj_list_size) {
      auto u = adj_list_v[index];
      if (ctx.dfs_timestamp[u] == SCCTarjanContext::NULL_INDEX) {
        // (b.1) implemented in a non-recursive manner
        ctx.dfs_stack.emplace(u);
        goto n_strongly_connected_components_dfs_non_recursive_head;
      } else if (ctx.in_trace_stack.test(u)) {
        // (b.3) Updates low[v] by dfn[u] on back edge detected
        ctx.low[v] = std::min(ctx.low[v], ctx.dfs_timestamp[u]);
      }
    }
    // (c) Counts for each DFS root
    if (ctx.low[v] == ctx.dfs_timestamp[v]) {
      ctx.result += 1;
      while (true) {
        BOOST_ASSERT(!ctx.trace_stack.empty());
        auto top = ctx.trace_stack.top();
        ctx.trace_stack.pop();
        ctx.in_trace_stack.reset(top);
        if (top == v) {
          break;
        }
      }
    }
    // (d) Returns
    BOOST_ASSERT(!ctx.dfs_stack.empty() && ctx.dfs_stack.top().v == v);
    ctx.dfs_stack.pop();
  }
}
} // namespace details

template <class... Attributes>
inline auto n_strongly_connected_components(const AdjacencyList<Attributes...>& graph) {
  auto ctx = details::SCCTarjanContext(graph::num_vertices(graph));
  for (auto v : vertices(graph)) {
    if (ctx.dfs_timestamp[v] == details::SCCTarjanContext::NULL_INDEX) {
      details::n_strongly_connected_components_dfs_non_recursive(ctx, graph, v);
    }
  }
  return ctx.result;
}

template <class... Attributes>
inline auto n_weakly_connected_components(const AdjacencyList<Attributes...>& graph,
                                          const InvAdjacencyList<Attributes...>& inv_graph,
                                          std::span<const vertex_id_t> skipped_list = {}) -> vertex_id_t {
  auto n = graph::num_vertices(graph);
  auto visited = DynamicBitset(n);
  auto queue = std::queue<vertex_id_t>{};

  for (auto v : skipped_list) {
    BOOST_ASSERT_MSG(v < n, "v is out of range [0, n).");
    visited.set(v); // Skips
  }

  auto push_to_queue = [&](vertex_id_t v) {
    queue.push(v);
    visited.set(v);
  };
  auto pop_from_queue = [&]() {
    auto res = queue.front();
    queue.pop();
    return res;
  };
  auto not_visited_filter = FILTER_VIEW(!visited.test(_1));

  auto res = vertex_id_t{0};
  for (auto v : vertices(graph)) {
    if (visited.test(v)) {
      continue;
    }
    // Starts BFS from v
    push_to_queue(v);
    while (!queue.empty()) {
      auto cur = pop_from_queue();
      for (auto out_neighbor : graph[cur] | views::keys | not_visited_filter) {
        push_to_queue(out_neighbor);
      }
      for (auto in_neighbor : inv_graph[cur] | views::keys | not_visited_filter) {
        push_to_queue(in_neighbor);
      }
    }
    res += 1; // New component {v, ...} counted
  }
  return res;
}

template <class... Attributes>
inline auto n_weakly_connected_components(const AdjacencyList<Attributes...>& graph,
                                          const InvAdjacencyList<Attributes...>& inv_graph,
                                          std::initializer_list<vertex_id_t> skipped_list) -> vertex_id_t {
  return n_weakly_connected_components(graph, inv_graph, std::span<const vertex_id_t>{skipped_list});
}
