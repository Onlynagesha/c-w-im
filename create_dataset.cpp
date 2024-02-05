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
auto get_p_values_generic(const ReadGraphParams& params, double in_deg, double out_deg,
                          double penalize_factor) noexcept {
  BOOST_ASSERT_MSG(in_deg > 0, "In-degree must be a positive integer.");
  BOOST_ASSERT_MSG(out_deg > 0, "Out-degree must be a positive integer.");
  BOOST_ASSERT_MSG(penalize_factor >= 1, "Penalization factor must fall in range [1, +inf)");

  auto _p0 = *params.lambda_in / in_deg / penalize_factor;
  auto _s0 = _p0 + *params.lambda_seed / out_deg / penalize_factor;
  auto _b0 = _p0 + *params.lambda_boost / out_deg / penalize_factor;

  auto clamp_fn = [&](auto p) { return std::clamp(p, *params.mu_min, *params.mu_max); };
  auto [p0, s0, b0] = tuple_transform(clamp_fn, std::tuple{_p0, _s0, _b0});

  edge_probability_t p = 1.0 - pow(1.0 - p0, *params.alpha);
  edge_probability_t p_seed = 1.0 - pow(1.0 - s0, *params.alpha_seed);
  edge_probability_t p_boost = 1.0 - pow(1.0 - b0, *params.alpha_boost);

#define CHECK_P_RANGE(z, _, p) BOOST_ASSERT_MSG(p >= 0 && p <= 1, BOOST_PP_STRINGIZE(p)" is out of range [0, 1]");
  BOOST_PP_SEQ_FOR_EACH(CHECK_P_RANGE, _, (p)(p_seed)(p_boost))

  if constexpr (std::is_same_v<E, WIMEdge>) {
    return WIMEdge{.p = p, .p_seed = p_seed};
  } else if constexpr (std::is_same_v<E, WBIMEdge>) {
    return WBIMEdge{.p = p, .p_seed = p_seed, .p_boost = p_boost};
  } else {
    static_assert(rfl::always_false_v<E>, "Invalid edge type.");
  }
}
} // namespace

auto ReadGraphParams::get_wim_p_values(double in_deg, double out_deg, double penalize_factor) const noexcept
    -> WIMEdge {
  return get_p_values_generic<WIMEdge>(*this, in_deg, out_deg, penalize_factor);
}

auto ReadGraphParams::get_wbim_p_values(double in_deg, double out_deg, double penalize_factor) const noexcept
    -> WBIMEdge {
  return get_p_values_generic<WBIMEdge>(*this, in_deg, out_deg, penalize_factor);
}

auto create_wim_dataset(const ReadGraphParams& params) noexcept -> WIMReadGraphResult {
  // 4 : indent as 4 spaces
  ELOGFMT(INFO, "Begins reading graph. Parameters: {:4}", params);
  auto raw_graph = graph::read_mm<DIRECTED>(params.input_file);
  graph::remove_self_loops(raw_graph);
  auto n = graph::num_vertices(raw_graph);
  auto m = graph::num_edges(raw_graph);

  auto out_degs = graph::degrees<0>(raw_graph);
  auto in_degs = graph::degrees<1>(raw_graph);
  auto sum_degree = LAMBDA_1(out_degs[_1] + in_degs[_1]);

  auto vertices_by_degree_all = [&] {
    auto view = range(n);
    return std::vector(view.begin(), view.end());
  }();
  ranges::sort(vertices_by_degree_all, ranges::less{}, sum_degree);

  auto vertex_weights = std::vector<vertex_weight_t>(n, 1.0_vw);
  // Only the vertices with less degrees are preserved
  auto n_weighted = static_cast<size_t>((1.0 - *params.unweighted_ratio) * n);
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
          make_histogram(vertex_weights, *params.histogram_width, *params.histogram_height));

  // Only the vertices with larger degrees are preserved
  auto n_penalized = static_cast<size_t>(*params.penalized_ratio * n);
  auto vertices_penalized = [&] {
    auto vertices = vertices_by_degree_all | views::drop(n - n_penalized);
    return boost::container::flat_set(vertices.begin(), vertices.end());
  }();
  if (n_penalized > 0 && *params.penalize_factor > 1.0) {
    ELOGFMT(INFO, "{} vertices ({:.3}% of all) with larger degree is penalized by {:.4}", n_penalized,
            *params.penalized_ratio * 100, *params.penalize_factor);
  }

  auto edges = DirectedEdgeList<WIMEdge>{};
  edges.open_for_push_back();
  for (auto [u, v] : raw_graph) {
    auto w = vertices_penalized.contains(u)
                 ? params.get_p_values<WIMEdge>(in_degs[v], out_degs[u], *params.penalize_factor)
                 : params.get_p_values<WIMEdge>(in_degs[v], out_degs[u]);
    edges.push_back(u, v, w);
  }
  edges.close_for_push_back();
  BOOST_ASSERT_MSG(m == graph::num_edges(edges), "# of edges mismatch.");

  ELOGFMT(INFO, "Finishes reading graph from file '{}'. |V| = {}, |E| = {}", params.input_file, n, m);
  return {.edge_list = std::move(edges), .vertex_weights = std::move(vertex_weights)};
}

auto create_wim_dataset_jstr(const std::string& params_json_text) noexcept -> rfl::Result<WIMReadGraphResult> {
  auto params = rfl::json::read<ReadGraphParams>(params_json_text);
  return params.transform(create_wim_dataset);
}
