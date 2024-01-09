#include "create_dataset.h"
#include "dump.h"
#include "utils/reflect.h"
#include <nlohmann/json.hpp>
#include <nwgraph/build.hpp>
#include <nwgraph/io/mmio.hpp>
#include <rfl.hpp>
#include <rfl/json.hpp>

#include <fmt/ranges.h>

auto create_wim_dataset(const ReadGraphParams& params) noexcept -> WIMReadGraphResult {
  // 4 : indent as 4 spaces
  fmt::println("Begins reading graph. Parameters: {:4}", params);
  auto raw_graph = graph::read_mm<DIRECTED>(params.input_file);
  graph::remove_self_loops(raw_graph);
  auto out_degs = graph::degrees<0>(raw_graph);
  auto in_degs = graph::degrees<1>(raw_graph);
  auto sum_degree = LAMBDA_1(out_degs[_1] + in_degs[_1]);

  auto edges = DirectedEdgeList<WIMEdge>{};
  edges.open_for_push_back();
  for (auto [u, v] : raw_graph) {
    auto w = params.get_p_values<WIMEdge>(in_degs[v], out_degs[u]);
    edges.push_back(u, v, w);
  }
  edges.close_for_push_back();

  auto vertex_weights = std::vector<vertex_weight_t>(graph::num_vertices(edges), 1.0_vw);
  auto vertices_by_degree = vertices(edges) | ranges::to<std::vector>();
  // Only the vertices with less degrees are preserved
  ranges::sort(vertices_by_degree, ranges::less{}, sum_degree);
  auto n_weighted = static_cast<ssize_t>(graph::num_vertices(edges) * params.weighted_ratio.value());
  vertices_by_degree.erase(vertices_by_degree.begin() + n_weighted, vertices_by_degree.end());
  // Let M = the maximum degree (in + out) among all the weighted vertices,
  // then for each vertex v, the expected weight of v is M / deg(v)
  if (!vertices_by_degree.empty()) {
    auto max_degree = ranges::max(vertices_by_degree | views::transform(sum_degree));
    for (auto v : vertices_by_degree | views::filter(LAMBDA_1(sum_degree(_1) < max_degree))) {
      auto d = sum_degree(v);
      // Decomposes to 1.0 + a geometric distribution whose expected value is M / deg(v) - 1
      auto dist = std::geometric_distribution<int64_t>{1.0 * d / max_degree};
      vertex_weights[v] = 1.0_vw + dist(rand_engine);
    }
  }

  constexpr auto msg_pattern = "Finishes reading graph from file '{}'. |V| = {}, |E| = {}";
  fmt::println(msg_pattern, params.input_file, graph::num_vertices(edges), graph::num_edges(edges));
  return {.edge_list = std::move(edges), .vertex_weights = std::move(vertex_weights)};
}

auto create_wim_dataset_jstr(const std::string& params_json_text) noexcept -> rfl::Result<WIMReadGraphResult> {
  auto params = rfl::json::read<ReadGraphParams>(params_json_text);
  return params.transform(create_wim_dataset);
}
