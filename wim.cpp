#include "wim.h"
#include "utils/easylog.h"
#include <fmt/ranges.h>

namespace {
template <same_as_either<double, std::vector<double>> ReturnType, //
          is_edge_property E,                                     //
          same_as_either<VertexSet, std::nullptr_t> BoostedSet,   //
          random_access_range_of_exactly<vertex_weight_t> VertexWeights>
auto wim_simulate_generic(const AdjacencyList<E>& graph, const VertexSet& seeds, const BoostedSet& boosted_vertices,
                          VertexWeights vertex_weights, uint64_t try_count) -> rfl::Result<ReturnType> {
  if (try_count <= 0) {
    return rfl::Error{"Try count of simulation must be a positive integer."};
  }
  auto n = graph::num_vertices(graph);
  auto s = seeds.vertex_list.size();
  auto sum = [&]() {
    if constexpr (std::is_same_v<ReturnType, double>) {
      return 0.0;
    } else {
      return std::vector<uint64_t>(n, 0);
    }
  }();

  // BFS queue: the first s vertices are fixed as seeds
  auto queue = make_reserved_vector<vertex_id_t>(n);
  ranges::copy(seeds.vertex_list, std::back_inserter(queue));
  auto vis = DynamicBitset(n);

  auto rand_test_edge = [&](vertex_id_t target, const E& e, bool source_is_seed) {
    if constexpr (std::is_same_v<E, WIMEdge>) {
      return e.rand_test(source_is_seed);
    } else if constexpr (std::is_same_v<E, WBIMEdge>) {
      if constexpr (std::is_same_v<BoostedSet, std::nullptr_t>) {
        return e.rand_test(false);
      } else {
        return e.rand_test(boosted_vertices.mask.test(target));
      }
    } else {
      static_assert(rfl::always_false_v<E>, "Unsupported edge property type.");
    }
  };
  auto add_to_queue = [&](vertex_id_t v) {
    vis[v] = true;
    queue.push_back(v);
  };

  for (auto try_index : range(try_count)) {
    // Initialization for current BFS process
    queue.resize(s); // Preserves seeds only
    vis = seeds.mask;
    // Starts BFS
    for (auto i = 0; i < queue.size(); i++) {
      auto cur = queue[i];
      for (auto [v, w] : graph[cur]) {
        if (!vis[v] && rand_test_edge(v, w, i < s)) {
          add_to_queue(v);
        }
      }
    }
    // Excluding the seeds themselves to prevent O(try_count * s) redundant calculation.
    if constexpr (std::is_same_v<ReturnType, double>) {
      auto non_seed_sum = accumulate_sum(queue | views::drop(s) | TRANSFORM_VIEW(vertex_weights[_1]));
      MYLOG_FMT_TRACE("Simulation #{}: score excluding seeds = {}, queue = {}", try_index, non_seed_sum, queue);
      sum += non_seed_sum;
    } else {
      ranges::for_each(queue | views::drop(s), LAMBDA_1(sum[_1] += 1));
      MYLOG_FMT_TRACE("After simulation #{}: sum = {}", try_index, sum);
    }
  }
  // Result = Total vertex weight of seeds + Average of total BFS-traversed vertex weight excluding seeds
  if constexpr (std::is_same_v<ReturnType, double>) {
    auto seed_sum = accumulate_sum(seeds.vertex_list | TRANSFORM_VIEW(vertex_weights[_1]));
    MYLOG_FMT_DEBUG("Total weight of seed vertices = {}, with # of seeds = {}", //
                    seed_sum, seeds.vertex_list.size());
    return accumulate_sum(seeds.vertex_list | TRANSFORM_VIEW(vertex_weights[_1])) + (sum / try_count);
  } else {
    auto res = [&]() {
      auto view = sum | TRANSFORM_VIEW(static_cast<double>(_1) / try_count);
      return std::vector(view.begin(), view.end());
    }();
    ranges::for_each(seeds.vertex_list, LAMBDA_1(res[_1] = 1.0));
    MYLOG_FMT_TRACE("Final result: {::.4f}", res);
    return res;
  }
}

template <is_edge_property E>
auto detect_probability_from_seeds_generic(const InvAdjacencyList<E>& inv_graph, const VertexSet& seeds,
                                           vertex_id_t max_distance) -> rfl::Result<DetectProbabilityFromSeedsResult> {
  auto n = graph::num_vertices(inv_graph);
  auto res = DetectProbabilityFromSeedsResult(n);
  auto temp = DetectProbabilityFromSeedsResult(n);

  auto paths = make_reserved_vector<edge_probability_t>(n);
  auto paths_boosted = make_reserved_vector<edge_probability_t>(n);
  for (auto k = 1_vid; k <= max_distance; k++) {
    for (auto v : vertices(inv_graph)) {
      if (seeds.contains(v)) {
        continue; // F[s] == F_boost[s] == 0.0 for each seed s
      }
      paths.clear();
      paths_boosted.clear();
      for (auto [u, w] : inv_graph[v]) {
        if (seeds.contains(u)) {
          paths.push_back(get_p_seed(w));
          paths_boosted.push_back(get_p_seed_or_boost(w));
        } else {
          paths.push_back(w.p * res.p_in[u]);
          paths_boosted.push_back(get_p_boost(w) * res.p_in[u]);
        }
      }
      temp.p_in[v] = at_least_1_probability_r(paths);
      temp.p_in_boosted[v] = at_least_1_probability_r(paths_boosted);
    }
    MYLOG_FMT_TRACE("F[{}][] = {::.4f}", k, temp.p_in);
    MYLOG_FMT_TRACE("F_boost[{}][] = {::.4f}", k, temp.p_in_boosted);
    if (temp.equals_with(res)) {
      MYLOG_FMT_TRACE("Early stops when k = {}", k);
      break;
    }
    temp.swap(res);
  }
  return std::move(res);
}
} // namespace

auto RRSketchSet::append_single(std::span<vertex_id_t> vertices) noexcept -> void {
  auto next_sketch = std::vector<vertex_id_t>{vertices.begin(), vertices.end()};
  auto next_sketch_id = sketches.size();
  for (auto v : vertices) {
    BOOST_ASSERT_MSG(v >= 0 && v < n_vertices(), "Vertex index out of range [0, n).");
    inv_sketches[v].push_back(next_sketch_id);
  }
  sketches.push_back(std::move(next_sketch));
}

auto RRSketchSet::append(size_t n_sketches) noexcept -> void {
  auto queue = make_reserved_vector<vertex_id_t>(n_vertices());
  auto queue_ext = make_reserved_vector<vertex_id_t>(n_vertices());
  auto vis = DynamicBitset{n_vertices()};
  auto vis_ext = DynamicBitset{n_vertices()};

  auto bfs_initialize = [](auto& queue, auto& vis, vertex_id_t center) {
    queue.assign({center});
    vis.reset();
    vis.set(center);
  };
  auto bfs_add_to_queue = [](auto& queue, auto& vis, vertex_id_t v) {
    BOOST_ASSERT_MSG(!vis.test(v), "Each vertex can not be added to BFS queue twice!");
    queue.push_back(v);
    vis.set(v);
  };

  for (auto _ : range(n_sketches)) {
    auto center = center_distribution(rand_engine);
    bfs_initialize(queue, vis, center);
    bfs_initialize(queue_ext, vis_ext, center);
    MYLOG_FMT_TRACE("Selects vertex #{} as the center of next RR-sketch.", center);

    for (auto i = 0zu; i < queue.size(); i++) {
      auto cur = queue[i];
      for (auto [v, w] : (*inv_graph)[cur]) {
        BOOST_ASSERT_MSG(w.is_valid(), "Invalid edge properties.");
        auto rand = rand_float();
        if (rand < w.p_seed && !vis_ext.test(v)) {
          // (1) with probability p_seed - p, v is added to RR-sketch yet unable to BFS further
          bfs_add_to_queue(queue_ext, vis_ext, v);
          // (2) with probability p, v is added to RR-sketch and BFS queue
          if (rand < w.p && !vis.test(v)) {
            bfs_add_to_queue(queue, vis, v);
          }
        }
      }
    }
    append_single(queue_ext);
  }
}

auto RRSketchSet::select(vertex_id_t k) noexcept -> std::vector<vertex_id_t> {
  k = std::min(k, n_vertices());
  // Uses negative count to mark the vertices that are already selected
  auto cover_count = [&] {
    auto view = inv_sketches | views::transform(ranges::ssize);
    return std::vector(view.begin(), view.end());
  }();
  auto covered = DynamicBitset{n_sketches()};
  auto res = std::vector<vertex_id_t>{};
  res.reserve(k);

  MYLOG_FMT_TRACE("Details of current RRSketchSet:\n{}", dump());

  for (auto i = 0_vid; i < k; i++) {
    MYLOG_FMT_TRACE("select: i = {}, cover_count = {}", i, cover_count);
    auto cur = static_cast<vertex_id_t>(ranges::max_element(cover_count) - cover_count.begin());
    res.push_back(cur);
    MYLOG_FMT_TRACE("Selects vertex #{} with cover_count = {}", cur, cover_count[cur]);
    BOOST_ASSERT_MSG(cover_count[cur] >= 0, "Invalid seed selected whose cover_count is negative!");

    cover_count[cur] = -1; // Marks as invalid
    for (auto r : inv_sketches[cur] | views::filter(LAMBDA_1(!covered.test(_1)))) {
      covered.set(r);
      for (auto v : sketches[r]) {
        cover_count[v] -= 1;
      }
    }
  }
  auto cover_percentage = 100.0 * covered.count() / n_sketches();
  MYLOG_FMT_DEBUG("Done selecting {} seed vertices. {:.3f}% of RR-sketch covered.", k, cover_percentage);
  return res;
}

auto RRSketchSet::dump() noexcept -> std::string {
  auto res = "RR-sketches:"s;
  for (auto [i, r] : views::enumerate(sketches)) {
    res += fmt::format("\n\t[{}] = {}", i, r);
  }
  res += "\nInverse of RR-sketches:";
  for (auto [v, ir] : views::enumerate(inv_sketches)) {
    res += fmt::format("\n\t[{}] = {}", v, ir);
  }
  return res;
}

auto wim_simulate(const AdjacencyList<WIMEdge>& graph, const VertexSet& seeds, uint64_t try_count) noexcept
    -> rfl::Result<double> {
  if (seeds.size() == 0) {
    return 0.0;
  }
  return wim_simulate_generic<double>(graph, seeds, nullptr, views::repeat(1.0_vw), try_count);
}

auto wim_simulate_w(const AdjacencyList<WIMEdge>& graph, std::span<const vertex_weight_t> vertex_weights,
                    const VertexSet& seeds, uint64_t try_count) noexcept -> rfl::Result<double> {
  if (vertex_weights.empty()) {
    return wim_simulate(graph, seeds, try_count);
  }
  return wim_simulate_generic<double>(graph, seeds, nullptr, vertex_weights, try_count);
}

auto wbim_simulate(const AdjacencyList<WBIMEdge>& graph, const VertexSet& seeds, const VertexSet& boosted_vertices,
                   uint64_t try_count) noexcept -> rfl::Result<double> {
  auto vertex_weights = views::repeat(1.0_vw);
  return wim_simulate_generic<double>(graph, seeds, boosted_vertices, vertex_weights, try_count);
}

auto wbim_simulate_w(const AdjacencyList<WBIMEdge>& graph, std::span<const vertex_weight_t> vertex_weights,
                     const VertexSet& seeds, const VertexSet& boosted_vertices, uint64_t try_count) noexcept
    -> rfl::Result<double> {
  if (vertex_weights.empty()) {
    return wbim_simulate(graph, seeds, boosted_vertices, try_count);
  }
  return wim_simulate_generic<double>(graph, seeds, boosted_vertices, vertex_weights, try_count);
}

auto wim_simulate_p(const AdjacencyList<WIMEdge>& graph, const VertexSet& seeds, uint64_t try_count)
    -> rfl::Result<std::vector<double>> {
  return wim_simulate_generic<std::vector<double>>(graph, seeds, nullptr, views::repeat(1.0_vw), try_count);
}

auto wbim_simulate_p(const AdjacencyList<WBIMEdge>& graph, const VertexSet& seeds, const VertexSet& boosted_vertices,
                     uint64_t try_count) -> rfl::Result<std::vector<double>> {
  return wim_simulate_generic<std::vector<double>>(graph, seeds, boosted_vertices, views::repeat(1.0_vw), try_count);
}

auto wim_detect_probability_from_seeds(const InvAdjacencyList<WIMEdge>& inv_graph, const VertexSet& seeds,
                                       vertex_id_t max_distance) -> rfl::Result<DetectProbabilityFromSeedsResult> {
  return detect_probability_from_seeds_generic(inv_graph, seeds, max_distance);
}

auto wbim_detect_probability_from_seeds(const InvAdjacencyList<WBIMEdge>& inv_graph, const VertexSet& seeds,
                                        vertex_id_t max_distance) -> rfl::Result<DetectProbabilityFromSeedsResult> {
  return detect_probability_from_seeds_generic(inv_graph, seeds, max_distance);
}
