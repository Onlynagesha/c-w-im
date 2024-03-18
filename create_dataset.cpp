#include "create_dataset.h"
#include "dump.h"
#include "utils/histogram.h"
#include "utils/reflect.h"
#include <boost/container/flat_set.hpp>
#include <boost/preprocessor/seq/for_each.hpp>
#include <boost/preprocessor/stringize.hpp>
#include <nlohmann/json.hpp>
#include <nwgraph/build.hpp>
#include <nwgraph/io/mmio.hpp>
#include <rfl.hpp>
#include <rfl/json.hpp>
#include <ylt/easylog.hpp>

namespace {
template <is_edge_property E>
auto create_dataset_generic(const typename CreateDatasetParamsTraits<E>::ParamsType& params) noexcept
    -> ReadGraphResult<E> {
  // 4 : indent as 4 spaces
  ELOGFMT(INFO, "Begins reading graph. Parameters: {:4}", params);
  auto raw_graph = graph::read_mm<DIRECTED>(params.common->input_file);
  graph::remove_self_loops(raw_graph);
  auto n = graph::num_vertices(raw_graph);
  auto m = graph::num_edges(raw_graph);

  auto out_degs = graph::degrees<0>(raw_graph);
  auto in_degs = graph::degrees<1>(raw_graph);
  auto sum_degree = LAMBDA_1(out_degs[_1] + in_degs[_1]);

  auto vertices_by_degree_all = std::vector<vertex_id_t>(n);
  std::iota(vertices_by_degree_all.begin(), vertices_by_degree_all.end(), 0);
  ranges::sort(vertices_by_degree_all, ranges::less{}, sum_degree);

  auto vertex_weights = std::vector<vertex_weight_t>(n, 1.0_vw);
  // Only the vertices with less degrees are preserved
  auto n_weighted = static_cast<size_t>((1.0 - *params.common->unweighted_ratio) * n);
  auto vertices_weighted = vertices_by_degree_all | views::take(n_weighted);
  // Let M = the maximum degree (in + out) among all the weighted vertices,
  // then for each vertex v, the expected weight of v is M / deg(v)
  if (!vertices_weighted.empty()) {
    auto max_degree = ranges::max(vertices_weighted | views::transform(sum_degree));
    for (auto v : vertices_weighted | views::filter(LAMBDA_1(sum_degree(_1) < max_degree))) {
      auto d = sum_degree(v);
      // Decomposes to 1.0 + a geometric distribution whose expected value is M / deg(v) - 1
      auto dist = std::geometric_distribution<int64_t>{1.0 * d / max_degree};
      vertex_weights[v] = 1.0_vw + dist(rand_engine);
    }
  }
  // Statistics of vertex weights
  ELOGFMT(INFO, "Distribution of vertex weights:\n{}",
          make_histogram(vertex_weights, *params.common->histogram_width, *params.common->histogram_height));

  // Only the vertices with larger degrees are preserved
  auto n_penalized = static_cast<size_t>(*params.common->penalized_ratio * n);
  auto vertices_penalized = [&] {
    auto vertices = vertices_by_degree_all | views::drop(n - n_penalized);
    return boost::container::flat_set(vertices.begin(), vertices.end());
  }();
  if (n_penalized > 0 && *params.common->penalize_factor > 1.0) {
    ELOGFMT(INFO, "{} vertices ({:.3}% of all) with larger degree is penalized by {:.4}", n_penalized,
            *params.common->penalized_ratio * 100, *params.common->penalize_factor);
  }

  auto edges = DirectedEdgeList<E>{};
  edges.open_for_push_back();
  for (auto [u, v] : raw_graph) {
    auto w = params.make_edge(in_degs[v], out_degs[u], vertices_penalized.contains(u));
    edges.push_back(u, v, w);
  }
  edges.close_for_push_back();
  BOOST_ASSERT_MSG(m == graph::num_edges(edges), "# of edges mismatch.");

  ELOGFMT(INFO, "Finishes reading graph from file '{}'. |V| = {}, |E| = {}", params.common->input_file, n, m);
  return {.edge_list = std::move(edges), .vertex_weights = std::move(vertex_weights)};
}
} // namespace

auto CreateWIMDatasetParams::make_edge(double in_deg, double out_deg, bool penalized) const noexcept -> WIMEdge {
  BOOST_ASSERT_MSG(in_deg > 0, "In-degree must be a positive integer.");
  BOOST_ASSERT_MSG(out_deg > 0, "Out-degree must be a positive integer.");

  auto penalize_factor_used = (penalized ? *common->penalize_factor : 1.0);
  auto _p0 = *common->lambda_in / in_deg / penalize_factor_used;
  auto _s0 = _p0 + *lambda_seed / out_deg / penalize_factor_used;

  auto clamp_fn = [&](auto p) { return std::clamp(p, *common->mu_min, *common->mu_max); };
  auto [p0, s0] = tuple_transform(clamp_fn, std::tuple{_p0, _s0});

  return WIMEdge{
      .p = static_cast<edge_probability_t>(1.0 - pow(1.0 - p0, *common->alpha)),
      .p_seed = static_cast<edge_probability_t>(1.0 - pow(1.0 - s0, *alpha_seed)),
  };
}

auto CreateWBIMDatasetParams::make_edge(double in_deg, double out_deg, bool penalized) const noexcept -> WBIMEdge {
  BOOST_ASSERT_MSG(in_deg > 0, "In-degree must be a positive integer.");
  BOOST_ASSERT_MSG(out_deg > 0, "Out-degree must be a positive integer.");

  auto penalize_factor_used = (penalized ? *common->penalize_factor : 1.0);
  auto _p0 = *common->lambda_in / in_deg / penalize_factor_used;
  auto _b0 = _p0 + *lambda_boost / out_deg / penalize_factor_used;

  auto clamp_fn = [&](auto p) { return std::clamp(p, *common->mu_min, *common->mu_max); };
  auto [p0, b0] = tuple_transform(clamp_fn, std::tuple{_p0, _b0});

  return WBIMEdge{
      .p = static_cast<edge_probability_t>(1.0 - pow(1.0 - p0, *common->alpha)),
      .p_boost = static_cast<edge_probability_t>(1.0 - pow(1.0 - b0, *alpha_boost)),
  };
}

auto create_wim_dataset(const CreateWIMDatasetParams& params) noexcept -> WIMReadGraphResult {
  return create_dataset_generic<WIMEdge>(params);
}

auto create_wim_dataset_jstr(const std::string& params_json_text) noexcept -> rfl::Result<WIMReadGraphResult> {
  return rfl::json::read<CreateWIMDatasetParams>(params_json_text).transform(create_wim_dataset);
}

auto create_wbim_dataset(const CreateWBIMDatasetParams& params) noexcept -> WBIMReadGraphResult {
  return create_dataset_generic<WBIMEdge>(params);
}

auto create_wbim_dataset_jstr(const std::string& params_json_text) noexcept -> rfl::Result<WBIMReadGraphResult> {
  return rfl::json::read<CreateWBIMDatasetParams>(params_json_text).transform(create_wbim_dataset);
}
