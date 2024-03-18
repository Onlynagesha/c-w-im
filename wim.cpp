#include "wim.h"
#include "dump.h"
#include "dump_utils.h"
#include "utils/easylog.h"
#include <deque>
#include <fmt/ranges.h>
#include <magic_enum_format.hpp>

namespace {
template <ranges::random_access_range SetRange, ranges::random_access_range InvSetRange>
  requires(forward_range_of<ranges::range_value_t<SetRange>, vertex_id_t> &&
           forward_range_of<ranges::range_value_t<InvSetRange>, size_t>)
auto greedy_max_cover(const SetRange& sets, const InvSetRange& inv_sets, vertex_id_t k) -> std::vector<vertex_id_t> {
  auto n_sets = ranges::size(sets);
  auto n_elements = static_cast<vertex_id_t>(ranges::size(inv_sets)); // Elements are indexed from 0 to n-1

  k = std::min(k, n_elements);
  // Uses negative count to mark the vertices that are already selected
  auto cover_count = [&] {
    auto view = inv_sets | views::transform(ranges::ssize);
    return std::vector(view.begin(), view.end());
  }();
  auto covered = DynamicBitset{n_sets};
  auto res = make_reserved_vector<vertex_id_t>(k);

  for (auto i = 0_vid; i < k; i++) {
    MYLOG_FMT_TRACE("greedy_max_cover before selecting result[{}]: cover_count = {}", i, cover_count);
    auto cur = static_cast<vertex_id_t>(ranges::max_element(cover_count) - cover_count.begin());
    res.push_back(cur);
    MYLOG_FMT_TRACE("Selects element #{} with cover_count = {}", cur, cover_count[cur]);
    BOOST_ASSERT_MSG(cover_count[cur] >= 0, "Invalid element selected whose cover_count is negative!");

    cover_count[cur] = -1; // Marks as invalid
    for (auto r : inv_sets[cur] | FILTER_VIEW(!covered.test(_1))) {
      covered.set(r);
      for (auto v : sets[r]) {
        cover_count[v] -= 1;
      }
    }
  }
  auto cover_percentage = 100.0 * covered.count() / n_sets;
  MYLOG_FMT_DEBUG("Done greedy-selecting {} elements. {:.3f}% of sets covered.", k, cover_percentage);
  return res;
}

template <same_as_either<double, std::vector<double>> ReturnType, //
          is_edge_property E,                                     //
          same_as_either<VertexSet, std::nullptr_t> BoostedSet,   //
          random_access_range_of_exactly<vertex_weight_t> VertexWeights>
auto simulate_generic(const AdjacencyList<E>& graph, const VertexSet& seeds, const BoostedSet& boosted_vertices,
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

auto RRSketchSet::select(vertex_id_t k) const noexcept -> std::vector<vertex_id_t> {
  return greedy_max_cover(sketches, inv_sketches, k);
}

auto RRSketchSet::dump() const noexcept -> std::string {
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

namespace {
struct PRRSketchQueueItem {
  vertex_id_t v;
  vertex_id_t dist;
};

struct PRRSketchEdge {
  vertex_id_t source;
  vertex_id_t target;
  PRRSketchEdgeState state;

  auto dump() const noexcept -> std::string {
    return fmt::format("{{.source = {}, .target = {}, .state = {}}}", source, target, state);
  }
};

struct PRRSketchCache {
  std::vector<vertex_id_t> dist_r;
  std::vector<vertex_id_t> dist_s;
  std::deque<PRRSketchQueueItem> queue;

  std::vector<vertex_id_t> seeds_found;
  std::vector<PRRSketchEdge> edges;

  std::vector<vertex_id_t> indices;
  DynamicBitset mask;
  DynamicBitset mask_dist_1;
  DynamicBitset mask_critical;

  auto dump(int indent = 0, int level = 0) const noexcept -> std::string {
    indent = std::max(indent, 0);
    level = std::max(level, 0);

    auto dynamic_bitset_to_string = [indent, level](const DynamicBitset& bitset) -> std::string {
      auto indices = range(bitset.size()) | FILTER_VIEW(bitset.test(_1));
      return dump_utils::dump_integer_array_as_braced_list(indices, indent, level + 1);
    };

    auto components = {
        fmt::format(".dist_r = {}", dump_utils::dump_integer_array_as_braced_list( //
                                        dist_r | views::transform(to_signed), indent, level + 1)),
        fmt::format(".dist_s = {}", dump_utils::dump_integer_array_as_braced_list( //
                                        dist_s | views::transform(to_signed), indent, level + 1)),
        fmt::format(".queue = {}", queue.empty() ? "(Empty std::deque)" : "(Non-empty std::deque)"),
        fmt::format(".seeds_found = {}", dump_utils::dump_integer_array_as_braced_list( //
                                             seeds_found, indent, level + 1)),
        fmt::format(".edges = {}", dump_utils::merge_dumped_components( //
                                       edges | views::transform(&PRRSketchEdge::dump), indent, level + 1)),
        fmt::format(".indices = {}", dump_utils::dump_integer_array_as_braced_list( //
                                         indices, indent, level + 1)),
        fmt::format(".mask = {}", dynamic_bitset_to_string(mask)),
        fmt::format(".mask_dist_1 = {}", dynamic_bitset_to_string(mask_dist_1)),
        fmt::format(".mask_critical = {}", dynamic_bitset_to_string(mask_critical)),
    };
    return dump_utils::merge_dumped_components(components, indent, level);
  }
};

auto prr_sketch_backward_bfs_1(const InvAdjacencyList<WBIMEdge>* inv_graph, const VertexSet* seeds, vertex_id_t center,
                               vertex_id_t k, PRRSketchCache* cache) -> bool {
  if (seeds->contains(center)) {
    MYLOG_FMT_TRACE("Skips PRR-sketch with center #{} which is a seed.", center);
    return false;
  }
  auto n = graph::num_vertices(*inv_graph);

  cache->dist_r.assign(n, PRRSketchSet::NULL_INDEX);
  cache->dist_r[center] = 0;
  cache->queue.assign({{.v = center, .dist = 0}});

  cache->edges.clear();
  cache->seeds_found.clear();

  while (!cache->queue.empty()) {
    [[maybe_unused]] auto queue_dist_diff = cache->queue.back().dist - cache->queue.front().dist;
    BOOST_ASSERT_MSG(queue_dist_diff == 0 || queue_dist_diff == 1,
                     "Unexpected behavior in BFS #1: monotonicity of queue violated.");

    auto [cur, cur_dist] = cache->queue.front();
    cache->queue.pop_front();
    if (cur_dist > cache->dist_r[cur]) {
      continue; // Visited before
    }
    BOOST_ASSERT_MSG(cur_dist < PRRSketchSet::NULL_INDEX, "Unexpected behavior: null index.");

    for (auto [v, w] : (*inv_graph)[cur]) {
      auto state = get_random_edge_state(w);
      if (state == PRRSketchEdgeState::BLOCKED) {
        continue; // Non-blocked edges accessible only
      }
      auto next_dist = cur_dist + (state == PRRSketchEdgeState::LIVE_UPON_BOOST);
      if (next_dist > k) {
        continue; // Pruning
      }
      cache->edges.push_back({.source = v, .target = cur, .state = state});
      if (next_dist >= cache->dist_r[v]) {
        continue;
      }
      cache->dist_r[v] = next_dist; // next_dist < dist[v]
      if (seeds->contains(v)) {
        cache->seeds_found.push_back(v);
        // Gain is always 0 if there's a LIVE path from some seed to the center
        if (next_dist == 0) {
          constexpr auto msg_pattern =
              "Skips PRR-sketch with center #{} since there's a LIVE path from seed #{} to center";
          MYLOG_FMT_TRACE(msg_pattern, center, v);
          return false;
        }
        // Otherwise, BFS stops here when reaching a seed vertex
      } else if (next_dist <= cur_dist) {
        cache->queue.push_front({.v = v, .dist = next_dist});
      } else {
        cache->queue.push_back({.v = v, .dist = next_dist});
      }
    }
  }
  // If not seed_fount, then there's no path from seed to center in current PRR-sketch
  if (cache->seeds_found.empty()) {
    MYLOG_FMT_TRACE("Skips PRR-sketch with center #{} since no seed vertex is reached.", center);
    return false;
  }
  MYLOG_FMT_TRACE("After backward BFS of step 1: cache = {:4}", *cache);
  return true; // OK
}

auto prr_sketch_build_uncompressed_sketch(vertex_id_t center, PRRSketchCache* cache) -> PRRSketch {
  cache->indices.clear();
  for (auto [u, v, w] : cache->edges) {
    cache->indices.push_back(u);
    cache->indices.push_back(v);
  }
  ranges::sort(cache->indices);
  auto [erase_first, erase_last] = ranges::unique(cache->indices);
  cache->indices.erase(erase_first, erase_last);

  auto n_of_sketch = static_cast<vertex_id_t>(cache->indices.size());
  auto vertices = FlatSet<vertex_id_t>(ORDERED_UNIQUE_RANGE, cache->indices.begin(), cache->indices.end());

  auto edge_list = DirectedEdgeList<PRRSketchEdgeState>{n_of_sketch};
  edge_list.open_for_push_back();
  // To mapped indices
  for (auto [u, v, w] : cache->edges) {
    BOOST_ASSERT(vertices.contains(u) && vertices.contains(v));
    u = static_cast<vertex_id_t>(vertices.find(u) - vertices.begin());
    v = static_cast<vertex_id_t>(vertices.find(v) - vertices.begin());
    edge_list.push_back(u, v, w);
  }
  edge_list.close_for_push_back();

  BOOST_ASSERT(vertices.contains(center));
  auto mapped_center = static_cast<vertex_id_t>(vertices.find(center) - vertices.begin());
  auto res = PRRSketch{
      .vertices = std::move(vertices),
      .critical_vertices = {}, // Leaves empty in the temporary uncompressed PRR-sketch
      .center = center,
      .mapped_center = mapped_center,
      .mapped_graph = {edge_list},
      .inv_mapped_graph = {edge_list},
  };
  MYLOG_FMT_TRACE("After building uncompressed PRR-sketch: cache = {:4}", *cache);
  MYLOG_FMT_TRACE("Uncompressed PRR-sketch = {:4}", res);
  return res;
}

// queue and dist should be initialized in prior
template <int IsInv, class Boosted = std::nullptr_t>
  requires(same_as_either<Boosted, std::nullptr_t, DynamicBitset>)
auto prr_sketch_bfs_2_common(const graph::adjacency<IsInv, PRRSketchEdgeState>& graph,
                             std::deque<PRRSketchQueueItem>& queue, std::span<vertex_id_t> dist,
                             vertex_id_t k = std::numeric_limits<vertex_id_t>::max(),
                             const Boosted& boosted_vertices = nullptr) -> void {
  // If boosted_vertices is present, and target is a boosted vertex, then the edge is considered as LIVE
  auto dist_increment = [&boosted_vertices](vertex_id_t target, PRRSketchEdgeState state) {
    auto is_lub = (state == PRRSketchEdgeState::LIVE_UPON_BOOST);
    if constexpr (std::is_same_v<Boosted, std::nullptr_t>) {
      return is_lub;
    } else {
      return is_lub && !boosted_vertices.test(target);
    }
  };

  while (!queue.empty()) {
    [[maybe_unused]] auto queue_dist_diff = queue.back().dist - queue.front().dist;
    BOOST_ASSERT_MSG(queue_dist_diff == 0 || queue_dist_diff == 1,
                     "Unexpected behavior in BFS #2: monotonicity of queue violated.");

    auto [cur, cur_dist] = queue.front();
    queue.pop_front();
    if (cur_dist > dist[cur]) {
      continue; // Visited before
    }
    // If boosted_vertices is present, then:
    // (1) If IsInv == 0, then the edge cur --> v is considered as LIVE if v is boosted;
    // (2) If IsInv == 1, then the edge v --> cur is considered as LIVE if cur is boosted
    for (auto [v, w] : graph[cur]) {
      auto next_dist = cur_dist + dist_increment(IsInv ? cur : v, w);
      if (next_dist > k || next_dist >= dist[v]) {
        continue;
      }
      dist[v] = next_dist; // next_dist < dist[v]
      if (next_dist <= cur_dist) {
        queue.push_front({.v = v, .dist = next_dist});
      } else {
        queue.push_back({.v = v, .dist = next_dist});
      }
    }
  }
}

// sketch is uncompressed
auto prr_sketch_forward_bfs_2(const PRRSketch& sketch, PRRSketchCache* cache) -> void {
  auto n_in_sketch = graph::num_vertices(sketch.mapped_graph);
  cache->dist_s.assign(n_in_sketch, PRRSketchSet::NULL_INDEX);
  cache->queue.clear();
  // Note: all the vertices, including the seeds, are stored in sketch.non_seed_vertices
  // in the temporarily uncompressed sketch.
  for (auto s : cache->seeds_found) {
    auto iter = sketch.vertices.find(s);
    BOOST_ASSERT(iter != sketch.vertices.end());
    s = static_cast<vertex_id_t>(iter - sketch.vertices.begin()); // To mapped index

    cache->dist_s[s] = 0;
    cache->queue.push_back({.v = s, .dist = 0});
  }
  // dist_s[v] = Mimimum distance from any seed vertex to v
  prr_sketch_bfs_2_common(sketch.mapped_graph, cache->queue, cache->dist_s);
  MYLOG_FMT_TRACE("After forward BFS of step 2: {:4}", *cache);
}

// sketch is uncompressed
auto prr_sketch_backward_bfs_2(const PRRSketch& sketch, PRRSketchCache* cache) -> void {
  auto n_in_sketch = graph::num_vertices(sketch.mapped_graph);
  cache->dist_r.assign(n_in_sketch, PRRSketchSet::NULL_INDEX);
  cache->dist_r[sketch.mapped_center] = 0;
  cache->queue.assign({{.v = sketch.mapped_center, .dist = 0}});

  // dist_r[v] = Mimimum distance from v to center
  prr_sketch_bfs_2_common(sketch.inv_mapped_graph, cache->queue, cache->dist_r);
  MYLOG_FMT_TRACE("After backward BFS of step 2: {:4}", *cache);
}

auto prr_sketch_build_compressed_sketch(const PRRSketch& uncompressed, vertex_id_t k, PRRSketchCache* cache)
    -> PRRSketch {
  auto n_in_uncompressed_sketch = graph::num_vertices(uncompressed.mapped_graph);
  cache->mask.clear();
  cache->mask.resize(n_in_uncompressed_sketch);
  cache->indices.clear();
  // Records all the mapped vertices that will be preserved in the final compressed PRR-sketch
  for (auto v : vertices(uncompressed.mapped_graph)) {
    if (to_signed(cache->dist_s[v]) <= 0) {
      continue; // To be merged as super seed, or inaccessible from seed vertices
    }
    if (cache->dist_r[v] >= PRRSketchSet::NULL_INDEX || cache->dist_s[v] + cache->dist_r[v] > k) {
      continue; // Inaccessible to center, or pruned
    }
    cache->mask.set(v);
    cache->indices.push_back(v);
  }
  BOOST_ASSERT(ranges::is_sorted(cache->indices));

  auto get_compressed_index = [&](vertex_id_t v) {
    auto it = ranges::find(cache->indices, v);
    BOOST_ASSERT(it != cache->indices.end());
    // 1 : Starts from 1 since 0 is reserved by the "super seed"
    return static_cast<vertex_id_t>(it - cache->indices.begin()) + 1;
  };
  auto mapped_center = get_compressed_index(uncompressed.mapped_center);

  cache->mask_critical.clear();
  cache->mask_critical.resize(n_in_uncompressed_sketch);

  // 1 : Plus one "super seed"
  auto edge_list = DirectedEdgeList<PRRSketchEdgeState>{cache->indices.size() + 1};
  edge_list.open_for_push_back();
  // Part 1: "Super seed" -> out-neighbors
  // mask_dist_1 records all the vertices v that there's an edge from "super seed" to v
  cache->mask_dist_1.clear();
  cache->mask_dist_1.resize(n_in_uncompressed_sketch);
  for (auto uc : vertices(uncompressed.mapped_graph)) {
    if (cache->dist_s[uc] != 0) {
      continue; // Not a member of "super seed"
    }
    for (auto [vc, wc] : uncompressed.mapped_graph[uc]) {
      if (!cache->mask.test(vc)) {
        continue;
      }
      BOOST_ASSERT_MSG(cache->dist_s[vc] == 1, "Unexpected behavior.");
      BOOST_ASSERT_MSG(wc == PRRSketchEdgeState::LIVE_UPON_BOOST, "Otherwise, dist_s[vc] will be 0.");
      cache->mask_dist_1.set(vc);
      // The only case that vc is a critical vertex
      if (cache->dist_r[vc] == 0) {
        cache->mask_critical.set(vc);
      }
    }
  }
  // 0 : Index of the "super seed"
  for (auto uc : vertices(uncompressed.mapped_graph) | FILTER_VIEW(cache->mask_dist_1.test(_1))) {
    edge_list.push_back(0, get_compressed_index(uc), PRRSketchEdgeState::LIVE_UPON_BOOST);
  }
  // Part 2: Edges whose source is not merged into the "super seed"
  for (auto uc : vertices(uncompressed.mapped_graph)) {
    if (!cache->mask.test(uc)) {
      continue; // Merged as "super seed", or inaccessible, or pruned
    }
    auto u = get_compressed_index(uc);
    // Links all the vertices whose dist_r == 0 to center directly, excluding the self loop
    if (cache->dist_r[uc] == 0) {
      if (uc != uncompressed.mapped_center) {
        edge_list.push_back(u, mapped_center, PRRSketchEdgeState::LIVE); // Links to center directly
      }
      continue;
    }
    for (auto [vc, w] : uncompressed.mapped_graph[uc]) {
      if (!cache->mask.test(vc)) {
        continue; // Merged as "super seed", or inaccessible, or pruned
      }
      auto v = get_compressed_index(vc);
      edge_list.push_back(u, v, w);
    }
  }
  edge_list.close_for_push_back();

  auto to_vertices_before_mapping = [&](const DynamicBitset& bitset) {
    auto head = uncompressed.vertices.begin();
    auto view = range(n_in_uncompressed_sketch) | FILTER_VIEW(bitset.test(_1)) | TRANSFORM_VIEW(head[_1]);
    BOOST_ASSERT(ranges::is_sorted(view));
    return FlatSet<vertex_id_t>(ORDERED_UNIQUE_RANGE, view.begin(), view.end());
  };
  auto res = PRRSketch{
      .vertices = to_vertices_before_mapping(cache->mask),
      .critical_vertices = to_vertices_before_mapping(cache->mask_critical),
      .center = uncompressed.center,
      .mapped_center = mapped_center,
      .mapped_graph = {edge_list},
      .inv_mapped_graph = {edge_list},
  };
  MYLOG_FMT_TRACE("Finally, cache = {:4}", *cache);
  MYLOG_FMT_TRACE("Compressed PRR-sketch = {:4}", res);
  return res;
}

auto prr_sketch_add_coverage(const PRRSketch& sketch, size_t sketch_index,
                             std::span<const vertex_id_t> boosted_vertices, std::span<std::vector<size_t>> coverage,
                             PRRSketchCache* cache) -> void {
  auto n_in_sketch = graph::num_vertices(sketch.mapped_graph);

  // Original index -> Mapped index
  auto& boosted_set = cache->mask; // Resource reusing
  boosted_set.clear();
  boosted_set.resize(n_in_sketch);
  for (auto b : boosted_vertices) {
    if (auto idx = sketch.mapped_vertex_index(b); idx != PRRSketch::NULL_INDEX) {
      boosted_set.set(idx);
    }
  }

  // (1) Forward BFS (constrained such that dist_s <= 1)
  cache->dist_s.assign(n_in_sketch, PRRSketchSet::NULL_INDEX);
  cache->dist_s[0] = 0;                       // 0 : Index of the "super seed"
  cache->queue.assign({{.v = 0, .dist = 0}}); // 0 : Index of the "super seed"
  prr_sketch_bfs_2_common(sketch.mapped_graph, cache->queue, cache->dist_s, 1, boosted_set);

  // (2) Backward BFS (constrained such that dist_r <= 0)
  cache->dist_r.assign(n_in_sketch, PRRSketchSet::NULL_INDEX);
  cache->dist_r[sketch.mapped_center] = 0;
  cache->queue.assign({{.v = sketch.mapped_center, .dist = 0}});
  prr_sketch_bfs_2_common(sketch.inv_mapped_graph, cache->queue, cache->dist_r, 0, boosted_set);

  cache->mask_critical.clear();
  cache->mask_critical.resize(n_in_sketch);
  for (auto u : range(n_in_sketch)) {
    if (cache->dist_s[u] != 0) {
      continue;
    }
    for (auto [v, w] : sketch.mapped_graph[u]) {
      if (cache->dist_s[v] == 1 && cache->dist_r[v] == 0) {
        BOOST_ASSERT_MSG(w == PRRSketchEdgeState::LIVE_UPON_BOOST, "Otherwise, dist_s[v] == 0");
        cache->mask_critical.set(v);
      }
    }
  }
  MYLOG_FMT_TRACE("Cache after picking critical vertices = {:4}", *cache);

  auto vertices_arr = sketch.vertices.begin();
  for (auto v : views::iota(1_vid, n_in_sketch) | FILTER_VIEW(cache->mask_critical.test(_1))) {
    auto v_unmapped = vertices_arr[v - 1];
    coverage[v_unmapped].push_back(sketch_index);
  }
}
} // namespace

auto PRRSketchSet::append_single(vertex_id_t center, vertex_id_t k, void* cache_raw_ptr) noexcept -> bool {
  auto* cache = static_cast<PRRSketchCache*>(cache_raw_ptr);
  if (!prr_sketch_backward_bfs_1(inv_graph, seeds, center, k, cache)) {
    return false;
  }
  auto uncompressed = prr_sketch_build_uncompressed_sketch(center, cache);
  prr_sketch_forward_bfs_2(uncompressed, cache);
  prr_sketch_backward_bfs_2(uncompressed, cache);
  sketches.push_back(prr_sketch_build_compressed_sketch(uncompressed, k, cache));
  // inv_critical
  auto cur_sketch_id = sketches.size() - 1;
  for (auto c : sketches.back().critical_vertices) {
    inv_critical[c].push_back(cur_sketch_id);
  }
  return true;
}

auto PRRSketchSet::append(size_t n_sketches, vertex_id_t k) noexcept -> void {
  auto cache_obj = PRRSketchCache{};
  for (auto success_count = 0zu; success_count < n_sketches;) {
    total_n_attempts += 1;
    auto center = center_distribution(rand_engine);
    if (append_single(center, k, &cache_obj)) {
      MYLOG_FMT_TRACE("PRR-Sketch #{}: Picks vertex #{} as center.", success_count, center);
      success_count += 1;
    }
  }
}

auto PRRSketchSet::select(vertex_id_t k) const noexcept -> std::vector<vertex_id_t> {
  if (k <= 0) {
    return {}; // Corner case: Empty list if k == 0
  }
  auto n = n_vertices();
  k = std::min(k, n);

  auto res = make_reserved_vector<vertex_id_t>(k);
  auto best_critical = ranges::max(range(n), ranges::less{}, LAMBDA_1(inv_critical[_1].size()));
  MYLOG_FMT_DEBUG("Selects result[0] = {} (which covers {} PRR-sketches)", //
                  best_critical, inv_critical[best_critical].size());
  res.push_back(best_critical);

  auto covered = DynamicBitset{n_sketches()};
  for (auto s : inv_critical[best_critical]) {
    covered.set(s);
  }
  auto coverage = std::vector<std::vector<size_t>>(n);

  constexpr auto dump_coverage = [](std::span<const std::vector<size_t>> coverage) {
    auto res = "{"s;
    for (const auto& [i, c] : coverage | views::enumerate) {
      res += fmt::format("\n\t[{}] = {}", i, c);
    }
    res += "\n}";
    return res;
  };

  auto cache = PRRSketchCache{};
  for (auto candidate_index = 1_vid; candidate_index < k; candidate_index++) {
    ranges::for_each(coverage, &std::vector<size_t>::clear);
    for (auto [sketch_index, sketch] : sketches | views::enumerate) {
      if (covered.test(sketch_index)) {
        continue;
      }
      prr_sketch_add_coverage(sketch, sketch_index, res, coverage, &cache);
      MYLOG_FMT_TRACE("After PRR-sketch #{}: coverage => {}", sketch_index, dump_coverage(coverage));
    }
    for (auto v : res) {
      coverage[v].assign({PRRSketchSet::NULL_INDEX}); // Marks as invalid
    }
    auto best = ranges::max(range(n), ranges::less{}, [&](vertex_id_t v) {
      if (!coverage[v].empty() && coverage[v].front() == PRRSketchSet::NULL_INDEX) {
        return ssize_t{-1}; // Selected before
      }
      return ranges::ssize(coverage[v]);
    });
    MYLOG_FMT_DEBUG("Selects result[{}] = {} (which covers {} more PRR-sketches)", //
                    candidate_index, best, coverage[best].size());
    res.push_back(best);
    for (auto sketch_index : coverage[best]) {
      covered.set(sketch_index);
    }
  }
  return res;
}

auto PRRSketchSet::select_by_critical(vertex_id_t k) const noexcept -> std::vector<vertex_id_t> {
  return greedy_max_cover(sketches | TRANSFORM_VIEW(std::span{_1.critical_vertices}), inv_critical, k);
}

auto wim_simulate(const AdjacencyList<WIMEdge>& graph, const VertexSet& seeds, uint64_t try_count) noexcept
    -> rfl::Result<double> {
  if (seeds.size() == 0) {
    return 0.0;
  }
  return simulate_generic<double>(graph, seeds, nullptr, views::repeat(1.0_vw), try_count);
}

auto wim_simulate_w(const AdjacencyList<WIMEdge>& graph, std::span<const vertex_weight_t> vertex_weights,
                    const VertexSet& seeds, uint64_t try_count) noexcept -> rfl::Result<double> {
  if (vertex_weights.empty()) {
    return wim_simulate(graph, seeds, try_count);
  }
  return simulate_generic<double>(graph, seeds, nullptr, vertex_weights, try_count);
}

auto wbim_simulate(const AdjacencyList<WBIMEdge>& graph, const VertexSet& seeds, const VertexSet& boosted_vertices,
                   uint64_t try_count) noexcept -> rfl::Result<double> {
  auto vertex_weights = views::repeat(1.0_vw);
  return simulate_generic<double>(graph, seeds, boosted_vertices, vertex_weights, try_count);
}

auto wbim_simulate_w(const AdjacencyList<WBIMEdge>& graph, std::span<const vertex_weight_t> vertex_weights,
                     const VertexSet& seeds, const VertexSet& boosted_vertices, uint64_t try_count) noexcept
    -> rfl::Result<double> {
  if (vertex_weights.empty()) {
    return wbim_simulate(graph, seeds, boosted_vertices, try_count);
  }
  return simulate_generic<double>(graph, seeds, boosted_vertices, vertex_weights, try_count);
}

auto wim_simulate_p(const AdjacencyList<WIMEdge>& graph, const VertexSet& seeds, uint64_t try_count)
    -> rfl::Result<std::vector<double>> {
  return simulate_generic<std::vector<double>>(graph, seeds, nullptr, views::repeat(1.0_vw), try_count);
}

auto wbim_simulate_p(const AdjacencyList<WBIMEdge>& graph, const VertexSet& seeds, const VertexSet& boosted_vertices,
                     uint64_t try_count) -> rfl::Result<std::vector<double>> {
  return simulate_generic<std::vector<double>>(graph, seeds, boosted_vertices, views::repeat(1.0_vw), try_count);
}

auto wbim_activation_probability_from_seeds(const InvAdjacencyList<WBIMEdge>& inv_graph, const VertexSet& seeds,
                                            vertex_id_t max_distance) -> rfl::Result<WBIMActivationProbability> {
  auto n = graph::num_vertices(inv_graph);
  auto res = WBIMActivationProbability(n);
  auto temp = WBIMActivationProbability(n);

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
          paths.push_back(w.p);
          paths_boosted.push_back(w.p_boost);
        } else {
          paths.push_back(w.p * res.p_in[u]);
          paths_boosted.push_back(w.p_boost * res.p_in[u]);
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
