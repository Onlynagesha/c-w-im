#include "wim.h"
#include "utils/easylog.h"
#include <fmt/ranges.h>

namespace {
template <is_edge_property E, same_as_either<VertexSet, std::nullptr_t> BoostedSet, class VertexWeights>
auto wim_simulate_generic(const AdjacencyList<E>& graph, const VertexSet& seeds, const BoostedSet& boosted_vertices,
                          VertexWeights vertex_weights, uint64_t try_count) -> rfl::Result<double> {
  if (try_count <= 0) {
    return rfl::Error{"Try count of simulation must be a positive integer."};
  }
  auto n = graph::num_vertices(graph);
  auto s = seeds.vertex_list.size();
  auto sum = 0.0;

  // BFS queue: the first s vertices are fixed as seeds
  auto queue = make_reserved_vector<vertex_id_t>(n);
  ranges::copy(seeds.vertex_list, std::back_inserter(queue));
  auto vis = DynamicBitset(n);

  auto rand_test_edge = [&](vertex_id_t source, const E& e, bool is_seed) {
    if constexpr (std::is_same_v<E, WIMEdge>) {
      return e.rand_test(is_seed);
    } else if constexpr (std::is_same_v<E, WBIMEdge>) {
      if constexpr (std::is_same_v<BoostedSet, std::nullptr_t>) {
        return e.rand_test(is_seed, false);
      } else {
        return e.rand_test(is_seed, boosted_vertices.mask.test(source));
      }
    } else {
      static_assert(rfl::always_false_v<E>, "Unsupported edge property type.");
    }
  };
  auto add_to_queue = [&](vertex_id_t v) {
    vis[v] = true;
    queue.push_back(v);
  };
  auto vid_to_weight = views::transform(LAMBDA_1(vertex_weights[_1]));

  for (auto try_index : range(try_count)) {
    // Initialization for current BFS process
    queue.resize(s); // Preserves seeds only
    vis = seeds.mask;
    // Part 1: Seed vertices (probability = p_seed)
    for (auto cur : seeds.vertex_list) {
      for (auto [v, w] : graph[cur]) {
        if (!vis[v] && rand_test_edge(cur, w, true)) {
          add_to_queue(v);
        }
      }
    }
    // Part 2: Non-seed vertices (probability = p or p_boost)
    for (auto i = s; i < queue.size(); i++) {
      auto cur = queue[i];
      for (auto [v, w] : graph[cur]) {
        if (!vis[v] && rand_test_edge(cur, w, false)) {
          add_to_queue(v);
        }
      }
    }
    // Excluding the seeds themselves to prevent O(try_count * s) redundant calculation.
    auto non_seed_sum = accumulate_sum(queue | views::drop(s) | vid_to_weight, 0.0);
    MYLOG_FMT_TRACE("Simulation #{}: score excluding seeds = {}, queue = {}", try_index, non_seed_sum, queue);
    sum += non_seed_sum;
  }
  // Result = Total vertex weight of seeds + Average of total BFS-traversed vertex weight excluding seeds
  auto seed_sum = accumulate_sum(seeds.vertex_list | vid_to_weight, 0.0);
  MYLOG_FMT_DEBUG("Total weight of seed vertices = {}, with # of seeds = {}", seed_sum, seeds.vertex_list.size());
  return accumulate_sum(seeds.vertex_list | vid_to_weight, 0.0) + (sum / try_count);
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
  return wim_simulate_generic(graph, seeds, nullptr, views::repeat(1.0_vw), try_count);
}

auto wim_simulate_w(const AdjacencyList<WIMEdge>& graph, std::span<const vertex_weight_t> vertex_weights,
                    const VertexSet& seeds, uint64_t try_count) noexcept -> rfl::Result<double> {
  return wim_simulate_generic(graph, seeds, nullptr, vertex_weights, try_count);
}

auto wbim_simulate(const AdjacencyList<WBIMEdge>& graph, const VertexSet& seeds, const VertexSet& boosted_vertices,
                   uint64_t try_count) noexcept -> rfl::Result<double> {
  auto vertex_weights = views::repeat(1.0_vw);
  auto before = wim_simulate_generic(graph, seeds, nullptr, vertex_weights, try_count);
  auto after = wim_simulate_generic(graph, seeds, boosted_vertices, vertex_weights, try_count);
  return after - before;
}

auto wbim_simulate_w(const AdjacencyList<WBIMEdge>& graph, std::span<const vertex_weight_t> vertex_weights,
                     const VertexSet& seeds, const VertexSet& boosted_vertices, uint64_t try_count) noexcept
    -> rfl::Result<double> {
  auto before = wim_simulate_generic(graph, seeds, nullptr, vertex_weights, try_count);
  auto after = wim_simulate_generic(graph, seeds, boosted_vertices, vertex_weights, try_count);
  return after - before;
}
