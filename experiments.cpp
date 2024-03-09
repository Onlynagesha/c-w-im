#include "experiments.h"
#include "dump.h"
#include "graph_connectivity.h"
#include "utils/easylog.h"
#include "utils/histogram.h"
#include "wim.h"
#include <fmt/ranges.h>
#include <magic_enum.hpp>
#include <nlohmann/json.hpp>

#define ASSERT_JSON_ARRAY_PTR(p) \
  BOOST_ASSERT_MSG((p) == nullptr || (p)->empty() || (p)->is_array(), "Expects " #p " to be an JSON array.")
#define ASSERT_JSON_OBJECT_PTR(p) \
  BOOST_ASSERT_MSG((p) == nullptr || (p)->empty() || (p)->is_object(), "Expects " #p " to be an JSON object.")

namespace {
using VertexList = std::vector<vertex_id_t>;
using VertexListList = std::vector<VertexList>;

auto init_easylog(const CommonExperimentParams& params) {
  easylog::set_min_severity(params.log_severity);
  // Note: async=true may trigger a bug that the program fails to terminate after everything is finished.
  if (!params.log_output_file.empty()) {
    easylog::init_log(params.log_severity, params.log_output_file, false, params.log_console);
  }
}

// Returns nullptr if out_path is empty
auto create_json_fout_ptr(const std::string& out_path) -> std::unique_ptr<std::ofstream> {
  auto json_fout_ptr = std::unique_ptr<std::ofstream>{};
  if (!out_path.empty()) {
    auto fout = std::ofstream{out_path};
    if (!fout.is_open()) {
      ELOGFMT(ERROR, "Failed to open JSON output file '{}'. Uses log output as fallback.", out_path);
    } else {
      ELOGFMT(INFO, "Successfully opens JSON output file '{}'.", out_path);
      json_fout_ptr = std::make_unique<std::ofstream>(std::move(fout));
    }
  } else {
    ELOG_INFO << "No JSON output file specified. Uses log output as fallback.";
  }
  return json_fout_ptr;
}

auto dump_to_json_fout_str(std::ofstream* fout, std::string_view contents) -> void {
  if (fout != nullptr) {
    MYLOG_FMT_DEBUG("JSON output: {}", contents);
    (*fout) << contents;
  } else {
    ELOGFMT(INFO, "JSON output: {}", contents);
  }
}

template <is_edge_property E>
auto dump_seeds_selected(const AdjacencyListPair<E>& graph, std::span<const vertex_id_t> seeds_selected)
    -> std::string {
  auto [n, m] = graph.graph_n_m();
  auto res = fmt::format("Seeds selected, with |E| / |V| = {:.3f}:", 1.0 * m / n);
  for (auto seed : seeds_selected) {
    BOOST_ASSERT_MSG(seed < n, "Seed index out of range [0, n).");
    res += fmt::format("\n\tvertex-id = {}, in-degree = {}, out-degree = {}", //
                       seed, graph.in_degree(seed), graph.out_degree(seed));
  }
  return res;
}

auto json_object_property(json* json_root, std::string_view key) -> json* {
  return (json_root != nullptr) ? &(*json_root)[key] : nullptr;
}

/*
json_root is expected to be an array.
(1) Accesses json_root[i]. If not existing, pads the array with empty objects until json_root[i] is created;
(2) Adds (or overwrites) all the items in kv_pairs_to_add to json_root[i];
(3) Returns the reference of json_root[i].

e.g. json_root = [ { “r": 1000, "value": 10 } ], kv_pairs_to_add = { "r": 10000, "s": 10 }
Then after the operation:
ith_json_array_element(json_root, 2, kv_pairs_to_add)["items"] = std::vector{1, 2, 3, 4};
json_root becomes: [
  { “r": 1000, "value": 10 },
  { }, // Padded empty object
  { "r": 10000, "s": 10, "items": [1, 2, 3, 4] }
] */
auto ith_json_array_element(json& json_root, size_t i, const json& kv_pairs_to_add) -> json& {
  BOOST_ASSERT_MSG(json_root.empty() || json_root.is_array(), "Non-array JSON object is not accepted.");
  // Fills with empty object until i == json_root.size() - 1
  auto& ith_element = [&]() -> json& {
    while (i >= json_root.size()) {
      json_root.emplace_back();
    }
    return json_root.back();
  }();
  for (const auto& [k, v] : kv_pairs_to_add.items()) {
    ith_element[k] = v; // May overwrite the existing properties
  }
  return ith_element;
}

auto ith_json_array_element(json* json_root, size_t i, const json& kv_pairs_to_add) -> json* {
  return (json_root != nullptr) ? &ith_json_array_element(*json_root, i, kv_pairs_to_add) : nullptr;
}

/*
json_root is expected to be an array.
(1) Accesses json_root[i]. If not existing, pads the array with empty objects until json_root[i] is created;
(2) Adds (or overwrites) all the items in kv_pairs_to_add to json_root[i];
(3) Accesses (or creates) json_root[i][key] and returns its reference.

e.g. json_root = [ { “r": 1000, "value": 10 } ], kv_pairs_to_add = { "r": 20000, "s": 10 }
Then after the operation:
ith_json_array_element_with_key(json_root, 2, "items", kv_pairs_to_add) = std::vector{1, 2, 3, 4};
json_root becomes: [
  { “r": 1000, "value": 10 },
  { }, // Padded empty object
  { "r": 20000, "s": 10, "items": [1, 2, 3, 4] }
] */
auto ith_json_array_element_with_key( //
    json& json_root, size_t i, std::string_view key, const json& kv_pairs_to_add) -> json& {
  // May create an empty object if not existing before
  return ith_json_array_element(json_root, i, kv_pairs_to_add)[key];
}

auto ith_json_array_element_with_key( //
    json* json_root, size_t i, std::string_view key, const json& kv_pairs_to_add) -> json* {
  return (json_root != nullptr) ? &ith_json_array_element_with_key(*json_root, i, key, kv_pairs_to_add) : nullptr;
}

/*
json_root will be added the following contents: {
  "n": 1500,    // # of vertices in the graph
  "m": 37500,   // # of directed edges in the graph
} */
template <is_edge_property E>
auto write_graph_basic_information(json& json_root, const AdjacencyList<E>& graph) {
  auto n = graph.num_vertices()[0];
  auto m = graph.num_edges();
  json_root["n"] = n;
  json_root["m"] = m;
  return std::tuple{n, m};
}

/*
json_root will be added the following contents: {
  "n_sketches": [10000, 100000],  // # of RR-sketches to be generated
  "n_seeds": [5, 10],             // # of seed vertices to be selected
  "simulation_try_count": 10000     // # of trials in monte-carlo simulation
} */
auto write_wim_params_information(json& json_root, const WIMParams& params) -> void {
  json_root["n_sketches"] = params.n_sketches;
  json_root["n_seeds"] = params.n_seeds;
  json_root["simulation_try_count"] = *params.simulation_try_count;
}

/*
json_root will be added the following contents: {
  "n_wcc": 12  // # of weakly connected components in the graph
} */
template <is_edge_property E>
auto write_graph_connectivity_information(json& json_root, const AdjacencyListPair<E>& graph) {
  // # of SCCs is absent since the recursive algorithm may fail due to DFS stack limit
  auto n_wcc = n_weakly_connected_components(graph.adj_list, graph.inv_adj_list);
  json_root["n_wcc"] = n_wcc;
  return n_wcc;
}

/*
Let res: SelectedSeedLists be the return value.
res[i] = List of length S. the seed vertices selected and sorted in descending order of estimated gain
when the # of RR-sketches is params.n_sketches[i],
where S = params.n_seeds.back() is the maximum number of seeds to be selected.

In the example below, all the 10 seeds selected are sorted by estimated gain in descending order.
The first 5 seeds, e.g. [0, 1275, 365, 45, 902], is the result when n_seeds = 5

Example of the resulted json_experiments (with n_seeds = 5, 10): [
  {
    “r": 1000,
    "average_sketch_size": 16.2,
    "selected_seeds": [ 0, 1275, 365, 45, 902, 125, 28, 692, 1234, 345 ]
  },
  {
    “r": 10000,
    "average_sketch_size": 16.5,
    "selected_seeds": [ 0, 902, 365, 1275, 45, 125, 1234, 28, 777, 345 ]
  }
]
Example of the resulted json_time_used: {
  "rr_sketch": [
    { “r": 1000, "seconds": 2.5 },  // Time usage of 10'000 RR-sketches in seconds
    { “r": 10000, "seconds": 24.0 } // ... 100'000 RR-sketches ...
  ],
  "select_seeds": [
    { “r": 1000, "seconds": 0.07 }, // Time usage of selecting 10 seeds with 10'000 RR-sketches
    { “r": 10000, "seconds": 0.56 } // ... with 100'000 RR-sketches
  ],
} */
auto do_wim_experiment_get_seeds(const WIMAdjacencyListPair& graph, const CommonExperimentParams& common_params,
                                 const WIMParams& wim_params, json* json_experiments = nullptr,
                                 json* json_time_used = nullptr) -> rfl::Result<VertexListList> {
  ASSERT_JSON_ARRAY_PTR(json_experiments);
  ASSERT_JSON_OBJECT_PTR(json_time_used);

  auto timer = nw::util::seconds_timer{};

  const auto& [adj_list, inv_adj_list, vertex_weights] = graph;
  auto [n, m] = graph_n_m(adj_list);
  auto max_n_seeds = ranges::max(wim_params.n_seeds);

  auto rr_sketches_total_time_used = 0.0;
  auto rr_sketches = RRSketchSet{&inv_adj_list, vertex_weights};

  auto* json_time_sketch = json_object_property(json_time_used, "rr_sketch");
  auto* json_time_select_seeds = json_object_property(json_time_used, "select_seeds");

  auto res = VertexListList{};
  res.reserve(wim_params.n_sketches.size());
  for (auto [i, r] : wim_params.n_sketches | views::enumerate) {
    auto r_new = r - rr_sketches.n_sketches();
    timer.start();
    rr_sketches.append(r_new);
    timer.stop();

    auto rr_sketches_new_time_used = timer.elapsed();
    rr_sketches_total_time_used += rr_sketches_new_time_used;

    auto avg_sketch_size = rr_sketches.average_sketch_size();
    // Log of time usage & statistics of RR-sketches
    constexpr auto msg_pattern_1 = //
        "Done appending new {1} RR-sketches in {3:.3f} sec."
        "\n\tAverage size of all the {0} RR-sketches = {4:.3f}"
        "\n\t\t(with {5:.2f}% of RR-sketches that contains 1 vertex only.)."
        "\n\tTotal time used = {2:.3f} sec."
        "\n\tEstimated memory usage = {6}.";
    ELOGFMT(INFO, msg_pattern_1,                                    //
            r, r_new,                                               // {0}, {1}
            rr_sketches_total_time_used, rr_sketches_new_time_used, // {2}, {3}
            avg_sketch_size,                                        // {4}
            rr_sketches.percentage_of_single_vertex_sketch(),       // {5}
            rr_sketches.rr_sketch_total_size_str()                  // {6}
    );
    // Histogram of the distribution of RR-sketch size
    MYLOG_FMT_DEBUG("Histogram of RR-sketch sizes:\n{}",
                    make_histogram(rr_sketches.sketch_sizes(), common_params.histogram_shape()));
    // JSON
    if (json_time_sketch != nullptr) {
      ith_json_array_element_with_key(*json_time_sketch, i, "seconds", {{"r", r}}) = rr_sketches_total_time_used;
    }
    if (json_experiments != nullptr) {
      ith_json_array_element_with_key(*json_experiments, i, "average_sketch_size", {{"r", r}}) = avg_sketch_size;
    }

    timer.start();
    auto seed_list = rr_sketches.select(max_n_seeds);
    timer.stop();
    ELOGFMT(INFO, "Done selecting {} seeds with {} RR-sketches. Time used = {:.3f} sec.", //
            max_n_seeds, r, timer.elapsed());
    MYLOG_DEBUG(dump_seeds_selected(graph, seed_list));

    // JSON output
    if (json_time_select_seeds != nullptr) {
      ith_json_array_element_with_key(*json_time_select_seeds, i, "seconds", {{"r", r}}) = timer.elapsed();
    }
    if (json_experiments != nullptr) {
      ith_json_array_element_with_key(*json_experiments, i, "seeds_selected", {{"r", r}}) = seed_list;
    }
    // Result
    res.push_back(std::move(seed_list));
  }
  return std::move(res);
}

// json_objective_function <- A floating-point value, F(S)
// json_time_used <- A floating-point value, time used in seconds
auto do_simulate(const WIMAdjacencyListPair& graph, std::span<const vertex_id_t> seeds_selected,
                 uint64_t simulation_try_count, json* json_objective_function = nullptr, json* json_time_used = nullptr)
    -> void {
  const auto& [adj_list, inv_adj_list, vertex_weights] = graph;
  auto n = graph::num_vertices(adj_list);
  auto timer = nw::util::seconds_timer{};
  timer.start();
  wim_simulate_w(adj_list, vertex_weights, {n, seeds_selected}, simulation_try_count)
      .and_then([&](double sim_res) -> ResultVoid {
        auto time_used = timer.lap();
        ELOGFMT(INFO, "Done simulation. Result = {:.3f}. Time used = {:.3f} sec.", sim_res, time_used);
        if (json_objective_function != nullptr) {
          *json_objective_function = sim_res;
        }
        if (json_time_used != nullptr) {
          *json_time_used = time_used;
        }
        return RESULT_VOID_SUCCESS;
      })
      .or_else([&](const rfl::Error& e) -> ResultVoid {
        // Reports the error to log output and continues the experiment
        ELOGFMT(ERROR, "\nFailed simulation with seed list {}.\n\tError message: `{}'", seeds_selected, e.what());
        return RESULT_VOID_SUCCESS;
      });
}

/*
Example of the resulted json_results: [
  { "s": 5, "objective_function": 632.4 },
  { "s": 10, "objective_function": 1172.5 }
]
Example of the resulted json_time_used: [
  { "s": 5, "seconds": 1.9 },
  { "s": 10, "seconds": 2.7 }
]
*/
auto do_simulate_with_n_seeds(const WIMAdjacencyListPair& graph, std::span<const vertex_id_t> seeds_selected,
                              std::span<const vertex_id_t> n_seeds, uint64_t simulation_try_count,
                              json* json_results = nullptr, json* json_time_used = nullptr) -> void {
  BOOST_ASSERT_MSG(ranges::max(n_seeds) <= ranges::size(seeds_selected), "No enough seeds!");
  ASSERT_JSON_ARRAY_PTR(json_results);
  ASSERT_JSON_ARRAY_PTR(json_time_used);

  for (auto [i, s] : n_seeds | views::enumerate) {
    auto* json_result_cur = ith_json_array_element_with_key(json_results, i, "objective_function", {{"s", s}});
    auto* json_time_cur = ith_json_array_element_with_key(json_time_used, i, "seconds", {{"s", s}});
    do_simulate(graph, seeds_selected.first(s), simulation_try_count, json_result_cur, json_time_cur);
  }
}

/*
Example of the resulted json_experiments: [
  { "r": 1000, "simulation_result": [ ... see above for details ... ] },
  { "r": 10000, "simulation_result": [ ... see above for details ... ] }
]
Example of the resulted json_time_used: [
  { "r": 1000, "by_n_seeds": [ ... see above for details ... ] },
  { "r": 10000, "by_n_seeds": [ ... see above for details ... ] }
] */
auto do_wim_experiment_simulate(const WIMAdjacencyListPair& graph, std::span<const VertexList> seeds_selected,
                                const WIMParams& params, json* json_experiments = nullptr,
                                json* json_time_used = nullptr) -> void {
  BOOST_ASSERT_MSG(seeds_selected.size() == params.n_sketches.size(), "Mismatch with # of RR-sketch groups.");
  ASSERT_JSON_ARRAY_PTR(json_experiments);
  ASSERT_JSON_ARRAY_PTR(json_time_used);

  for (auto [i, r] : params.n_sketches | views::enumerate) {
    auto* json_exp_cur = ith_json_array_element_with_key(json_experiments, i, "simulation_result", {{"r", r}});
    auto* json_time_cur = ith_json_array_element_with_key(json_time_used, i, "by_n_seeds", {{"r", r}});

    for (auto [j, s] : params.n_seeds | views::enumerate) {
      do_simulate(graph, std::span{seeds_selected[i]}.first(s), *params.simulation_try_count,
                  ith_json_array_element_with_key(json_exp_cur, j, "objective_function", {{"s", s}}),
                  ith_json_array_element_with_key(json_time_cur, j, "seconds", {{"s", s}}));
    }
  }
}

/*
Returns the coarsening brief info of each layer (which is used during seed expansion).

Example of the resulted json_root: {
  "coarsening_level": 5,  // # of coarsening levels
  "n_coarsened": 1500,    // # of vertices after coarsening
  "m_coarsened": 37500    // # of edges after coarsening
}
*json_time_used = (a floating point value, time used in seconds)
*/
auto do_continue_coarsening_graph(const WIMAdjacencyListPair& original_graph,
                                  std::vector<CoarsenGraphBriefResult>& coarsen_results, vertex_id_t destination_level,
                                  vertex_id_t coarsening_threshold, const CoarseningParams& params,
                                  json* json_root = nullptr, json* json_time_used = nullptr) -> void {
  ASSERT_JSON_OBJECT_PTR(json_root);

  auto initial_level = static_cast<vertex_id_t>(coarsen_results.size());
  if (initial_level >= destination_level) {
    return; // No need to coarsen the graph
  }
  auto n_coarsened = (initial_level == 0) ? original_graph.n_vertices() : coarsen_results.back().coarsened.n_vertices();

  auto timer = nw::util::seconds_timer{};
  timer.start();
  for (auto level = initial_level + 1; n_coarsened > coarsening_threshold && level <= destination_level; level++) {
    auto cur_result = [&]() {
      auto item = (level <= 1) ? coarsen_wim_graph_p(original_graph, params)
                               : coarsen_wim_graph_p(coarsen_results.back().coarsened, params);
      return &coarsen_results.emplace_back(std::move(item));
    }();
    n_coarsened = cur_result->details.n_coarsened;
    ELOGFMT(INFO, "Finishes coarsening level {}: |V|, |E| = {}", level, cur_result->coarsened.graph_n_m());
  }
  auto time_used = timer.lap();
  ELOGFMT(INFO, "Finishes coarsening from level {} to {}: Time used = {:.3f} sec.", //
          initial_level, destination_level, time_used);

  auto [nc, mc] = coarsen_results.back().coarsened.graph_n_m();
  if (json_root != nullptr) {
    (*json_root)["coarsening_level"] = destination_level;
    (*json_root)["n_coarsened"] = nc;
    (*json_root)["m_coarsened"] = mc;
  }
  if (json_time_used != nullptr) {
    *json_time_used = time_used;
  }
}

auto do_coarsen_graph(const WIMAdjacencyListPair& graph, vertex_id_t destination_level,
                      vertex_id_t coarsening_threshold, const CoarseningParams& params, json* json_root = nullptr,
                      json* json_time_used = nullptr) -> std::vector<CoarsenGraphBriefResult> {
  auto res = make_reserved_vector<CoarsenGraphBriefResult>(destination_level);
  do_continue_coarsening_graph(graph, res, destination_level, coarsening_threshold, params, json_root, json_time_used);
  return res;
}

/*
Returns the expanded seed list.

*json_expanded_seeds <- List of expanded seed vertices
*json_time_used <- A floating point value, time used in seconds
*/
auto do_expand_seeds(const WIMAdjacencyListPair& original_graph,
                     std::span<const CoarsenGraphBriefResult> coarsen_results, VertexList seeds,
                     vertex_id_t n_fast_expanding_levels, const ExpandingParams& raw_params,
                     json* json_expanded_seeds = nullptr, json* json_time_used = nullptr) -> rfl::Result<VertexList> {
  auto level = static_cast<vertex_id_t>(coarsen_results.size());
  auto rule_local = ExpandingParams{.seed_expanding_rule = SeedExpandingRule::LOCAL};
  const auto& [adj_list, inv_adj_list, vertex_weights] = original_graph;

  auto timer = nw::util::seconds_timer{};
  timer.start();
  // Rewinds to level 0 (i,e. seeds of the original graph)
  for (auto back_level : range(level) | views::reverse) {
    // Uses LOCAL policy for the top f levels, where f = n_fast_expanding_levels
    auto params_to_use = [&]() -> const ExpandingParams* {
      if (back_level < n_fast_expanding_levels) {
        ELOGFMT(INFO, "Forces to use LOCAL policy when expanding back to level {}.", back_level);
        return &rule_local;
      }
      ELOGFMT(INFO, "Uses {} policy when expanding back to level {}.", //
              magic_enum::enum_name(raw_params.seed_expanding_rule), back_level);
      return &raw_params;
    }();
    auto expand_res = [&]() {
      if (back_level == 0) {
        return expand_wim_seed_vertices_w(adj_list, vertex_weights, coarsen_results[0], seeds, *params_to_use);
      }
      auto& prev = coarsen_results[back_level - 1];
      auto& cur = coarsen_results[back_level];
      return further_expand_wim_seed_vertices(prev, cur, seeds, *params_to_use);
    }();
    RFL_RETURN_ON_ERROR(expand_res);
    MYLOG_FMT_DEBUG("Seeds expansion back to level {} -> {}", back_level, expand_res->expanded_seeds);
    seeds = expand_res->expanded_seeds; // Updates seeds by the expanded
  }
  auto time_used = timer.lap();

  if (json_expanded_seeds != nullptr) {
    *json_expanded_seeds = seeds;
  }
  if (json_time_used != nullptr) {
    *json_time_used = time_used;
  }
  return std::move(seeds);
}

/*
Returns the expanded seed lists.
Example of the resulted json_results: [
  { “r": 1000, "expanded_seeds": [ 0, 1275, 45, 365, 902, 28, 125, 1234, 345, 789 ] },
  { “r": 10000, "expanded_seeds": [ 0, 45, 902, 1275, 125, 789, 345, 28, 1504, 277 ] }
]
Example of the resulted json_time_used: "expand_seeds": [
  { “r": 1000, "seconds": 0.87 },
  { “r": 10000, "seconds": 1.02}
], */
auto do_coarsening_experiment_expand_seeds(const WIMAdjacencyListPair& original_graph,
                                           std::span<const CoarsenGraphBriefResult> coarsen_results,
                                           VertexListList seed_lists, std::span<const size_t> n_rr_sketches,
                                           vertex_id_t n_fast_expanding_levels, const ExpandingParams& raw_params,
                                           json* json_results = nullptr, json* json_time_used = nullptr)
    -> rfl::Result<VertexListList> {
  BOOST_ASSERT_MSG(seed_lists.size() == n_rr_sketches.size(), "Mismatch between # of experiment groups.");

  for (auto [i, cur_seed_list] : seed_lists | views::enumerate) {
    auto n_sketches = n_rr_sketches[i];
    MYLOG_FMT_DEBUG("Seed list #{}: Initial seeds with r = {}: {}", i, n_sketches, cur_seed_list);

    auto* json_expanded_seeds = ith_json_array_element_with_key(json_results, i, "expanded_seeds", {{"r", n_sketches}});
    auto* json_seconds = ith_json_array_element_with_key(json_time_used, i, "seconds", {{"r", n_sketches}});
    auto cur_expanded = do_expand_seeds(original_graph, coarsen_results, cur_seed_list, n_fast_expanding_levels,
                                        raw_params, json_expanded_seeds, json_seconds);
    RFL_RETURN_ON_ERROR(cur_expanded);
    cur_seed_list = *cur_expanded;
  }
  return std::move(seed_lists);
}

/*
Example of JSON output: {
  "n": 1500, "m": 37500, "n_sketches": [10000, 100000], "n_seeds": [5, 10], "simulation_try_count": 10000,
  "experiments": [
    { “r": 1000, "average_sketch_size": 16.2,
      "selected_seeds": [ 0, 1275, 365, 45, 902, 125, 28, 692, 1234, 345 ],
      "simulation_result": [
        { "s": 5, "objective_function": 632.4 },
        { "s": 10, "objective_function": 1172.5 },
      ]
    },
    { “r": 10000, ... }
  ],
  "time_used": {
    "rr_sketch": [ { “r": 1000, "seconds": 2.5 }, { “r": 10000, "seconds": 24.0 } ],
    "select_seeds": [ { “r": 1000, "seconds": 0.07 }, { “r": 10000, "seconds": 0.56 } ],
    "simulate": [
      { “r": 1000, "by_n_seeds": [ { "s": 5, "seconds": 1.6 }, { "s": 10, "seconds": 2.2 } ] },
      { “r": 10000, "by_n_seeds": [ { "s": 5, "seconds": 1.9 }, { "s": 10, "seconds": 2.7 } ] }
    ]
  }
} */
auto do_wim_experiment(const AdjacencyListPair<WIMEdge>& graph, const WIMExperimentParams& params)
    -> rfl::Result<json> {
  auto json_root = json{};
  auto& json_experiments = json_root["experiments"];
  auto& json_time_used = json_root["time_used"];

  write_graph_basic_information(json_root, graph.adj_list);
  write_wim_params_information(json_root, *params.wim);

  return do_wim_experiment_get_seeds(graph, *params.common, *params.wim, &json_experiments, &json_time_used)
      .transform([&](VertexListList seeds_selected) {
        do_wim_experiment_simulate(graph, seeds_selected, *params.wim, //
                                   &json_experiments["simulation_result"], &json_time_used["simulate"]);
        return std::move(json_root);
      });
}

/*
Example of JSON output: {
  "n": 1500, "m": 37500, "n_sketches": [10000, 100000], "n_seeds": [5, 10], "simulation_try_count": 10000,
  "experiments": [
    { "level": 0, "n": 750, "m": 21600, "rr_sketch": [
      { “r": 1000, "average_sketch_size": 16.2,
        "selected_seeds": [ 0, 1275, 365, 45, 902, 125, 28, 692, 1234, 345 ]
      },
      { “r": 10000, "average_sketch_size": 16.5, "selected_seeds": [ ... ] }
    ], "simulation": [
      { “r": 1000, "simulation_result": [
        { "s": 5, "objective_function": 625.4 },
        { "s": 10, "objective_function": 1142.5 },
      ] },
      { “r": 10000, "simulation_result": [
        { "s": 5, "objective_function": 675.3 },
        { "s": 10, "objective_function": 1316.4 },
      ] },
    ] },
    { "level": 1, "n_coarsened": 375, "m_coarsened": 14400, "rr_sketch": [
      { “r": 1000, "average_sketch_size": ..., "selected_seeds": [ ... ] },
      { “r": 10000, "average_sketch_size": ..., "selected_seeds": [ ... ] }
    ], "expanding_seeds": [ ... see above for details ], "simulation": [
      { “r": 1000, "simulation_result": [
        { "s": 5, "objective_function": 623.7 },
        { "s": 10, "objective_function": 1123.9 },
      ] },
      { “r": 10000, "simulation_result": [
        { "s": 5, "objective_function": 670.6 },
        { "s": 10, "objective_function": 1297.1 },
      ] },
    ] },
    { "level": 2, ... }
  ],
  "time_used": [
    { "level": 0,
      "rr_sketch": [ { “r": 1000, "seconds": 2.5 }, { “r": 10000, "seconds": 24.0 } ],
      "select_seeds": [ { “r": 1000, "seconds": 0.07 }, { “r": 10000, "seconds": 0.56 } ],
      "simulate": [
        { “r": 1000, "by_n_seeds": [ { "s": 5, "seconds": 1.6 }, { "s": 10, "seconds": 2.2 } ] },
        { “r": 10000, "by_n_seeds": [ { "s": 5, "seconds": 1.9 }, { "s": 10, "seconds": 2.7 } ] }
    },
    { "level": 1, "coarsen": 0.23,
      "rr_sketch": [ { “r": 1000, "seconds": ... }, { “r": 10000, "seconds": ... } ],
      "select_seeds": [ { “r": 1000, "seconds": ... }, { “r": 10000, "seconds": ... } ],
      "expand_seeds": [ ... see above for details ],
      "simulate": [
        { “r": 1000, "by_n_seeds": [ { "s": 5, "seconds": ... }, { "s": 10, "seconds": ... } ] },
        { “r": 10000, "by_n_seeds": [ { "s": 5, "seconds": ... }, { "s": 10, "seconds": ... } ] }
    },
  ]
} */
auto do_wim_coarsening_experiment(const AdjacencyListPair<WIMEdge>& graph, const WIMCoarseningExperimentParams& params)
    -> rfl::Result<json> {
  auto timer = nw::util::seconds_timer{};

  auto json_root = json{};
  auto& json_experiments_root = json_root["experiments"]; // As a list
  auto& json_time_used_root = json_root["time_used"];

  auto make_json_cur_items = [&](vertex_id_t level) -> std::tuple<json*, json*> {
    auto kv_pairs_to_add = json{{"level", level}};
    return std::tuple{&ith_json_array_element(json_experiments_root, level, kv_pairs_to_add),
                      &ith_json_array_element(json_time_used_root, level, kv_pairs_to_add)};
  };

  const auto& [adj_list, inv_adj_list, vertex_weights] = graph;
  write_wim_params_information(json_root, *params.wim);
  auto [n, m] = write_graph_basic_information(json_root, adj_list);

  if (n <= params.coarsening_threshold) {
    constexpr auto msg_pattern = "Experiments fails to proceed since |V| <= threshold, with |V| = {}, threshold = {}";
    return rfl::Error{fmt::format(msg_pattern, n, *params.coarsening_threshold)};
  }
  // Step 1: Level 0, i.e. initial graph
  ELOG_INFO << "Starts Level 0: solves WIM problem on the original graph.";
  {
    auto [json_exp_cur, json_time_cur] = make_json_cur_items(0);
    auto n_wcc = write_graph_connectivity_information(*json_exp_cur, graph);
    ELOGFMT(INFO, "# of weakly connected components = {}", n_wcc);

    RFL_RETURN_ON_ERROR(
        do_wim_experiment_get_seeds(graph, *params.common, *params.wim, &(*json_exp_cur)["rr_sketch"], json_time_cur)
            .transform([&](VertexListList seeds) {
              do_wim_experiment_simulate(graph, seeds, *params.wim, //
                                         &(*json_exp_cur)["simulation"], &(*json_time_cur)["simulate"]);
              return RESULT_VOID_SUCCESS;
            }));
  }
  // Step 2: Coarsening until the size threshold
  auto coarsen_results = std::vector<CoarsenGraphBriefResult>{};
  for (auto level = 1_vid, n_coarsened = n; n_coarsened > params.coarsening_threshold; ++level) {
    auto [json_exp_cur, json_time_cur] = make_json_cur_items(level);
    // Step 2.1: Gets the coarsened graph
    do_continue_coarsening_graph(graph, coarsen_results, level, *params.coarsening_threshold, *params.coarsening,
                                 json_exp_cur, &(*json_time_cur)["coarsen"]);
    const auto& cur_coarsen_result = coarsen_results.back();
    n_coarsened = cur_coarsen_result.coarsened.n_vertices();

    // Checks connectivity for debugging & data analysis
    auto n_wcc = write_graph_connectivity_information(*json_exp_cur, cur_coarsen_result.coarsened);
    ELOGFMT(INFO, "# of weakly connected components = {}", n_wcc);

    RFL_RETURN_ON_ERROR(
        // Step 2.2: Seed selection with the coarsened graph
        do_wim_experiment_get_seeds(cur_coarsen_result.coarsened, *params.common, *params.wim, //
                                    &(*json_exp_cur)["rr_sketch"], json_time_cur)
            .and_then([&](VertexListList seed_lists) -> rfl::Result<VertexListList> {
              // Step 2.3: Seed expansion to the original graph
              return do_coarsening_experiment_expand_seeds( //
                  graph, coarsen_results, std::move(seed_lists), params.wim->n_sketches,
                  *params.n_fast_expanding_levels, *params.expanding, &(*json_exp_cur)["expanding_seeds"],
                  &(*json_time_cur)["expand_seeds"]);
            })
            .transform([&](VertexListList expanded_seeds) {
              // Step 2.4: Simulation with expanded seeds
              do_wim_experiment_simulate(graph, expanded_seeds, *params.wim, //
                                         &(*json_exp_cur)["simulation"], &(*json_time_cur)["simulate"]);
              return RESULT_VOID_SUCCESS;
            }));
  } // for level = 1, 2, 3 ...
  return std::move(json_root);
}

auto do_wim_contrast_experiment(const AdjacencyListPair<WIMEdge>& graph, const WIMContrastExperimentParams& params)
    -> rfl::Result<json> {
  auto timer = nw::util::seconds_timer{};

  auto json_root = json{};
  auto& json_time_used = json_root["time_used"];

  json_root["n_seeds"] = params.n_seeds;
  json_root["simulation_try_count"] = *params.simulation_try_count;
  auto [n, m] = write_graph_basic_information(json_root, graph.adj_list);

  auto coarsen_results = do_coarsen_graph(graph, params.coarsening_level, params.coarsening_threshold,
                                          *params.coarsening, &json_root, &json_time_used["coarsen"]);
  auto* graph_for_algorithm = coarsen_results.empty() ? &graph : &coarsen_results.back().coarsened;
  auto n_experiments = 0;
  auto total_expanding_time = 0.0;

#define CONTRAST_EXPERIMENT_FN [&](const WIMAdjacencyListPair& graph) -> rfl::Result<VertexList>
  auto contrast_experiment_framework = [&](std::string_view algorithm_name, auto&& algorithm_fn) -> ResultVoid {
    ELOGFMT(INFO, "==== Starts contrast algorithm: {} ====", algorithm_name);
    auto& json_exp_cur = json_root[algorithm_name];
    auto& json_time_cur = json_time_used[algorithm_name];
    timer.start();
    // Step 1: Performs the algorithm
    return std::invoke(algorithm_fn, *graph_for_algorithm)
        .and_then([&](VertexList seeds) -> rfl::Result<VertexList> {
          timer.stop();
          ELOGFMT(INFO, "Finished algorithm '{}' in {:.3f} sec.", algorithm_name, timer.elapsed());
          MYLOG_FMT_DEBUG("Seed vertices before expansion = {}", seeds);
          json_time_cur["algorithm"] = timer.elapsed();

          // Step 2: Expands seed vertices
          return do_expand_seeds(graph, coarsen_results, seeds, *params.n_fast_expanding_levels, *params.expanding,
                                 &json_exp_cur["expanded_seeds"], &json_time_cur["expand"]);
        })
        .and_then([&](VertexList expanded_seeds) -> ResultVoid {
          // Step 3: Simulation with expanded seeds
          do_simulate_with_n_seeds(graph, expanded_seeds, params.n_seeds, *params.simulation_try_count, //
                                   &json_exp_cur["simulation"], &json_time_cur["simulate"]);
          return RESULT_VOID_SUCCESS;
        });
  };
  auto max_n_seeds = ranges::max(params.n_seeds);
  auto pagerank_params_common = PagerankParams{
      .damping_factor = *params.pagerank_damping_factor,
      .epsilon = *params.pagerank_epsilon,
      .k = max_n_seeds,
      .n_iterations = params.pagerank_n_iterations,
      .uses_vertex_weight = true,
      .uses_edge_weight = true,
  };
  if (params.with_max_degree) {
    RFL_RETURN_ON_ERROR(contrast_experiment_framework(
        "max_degree", CONTRAST_EXPERIMENT_FN { return max_out_degree(graph.adj_list, max_n_seeds); }));
  }
  if (params.with_max_strength) {
    RFL_RETURN_ON_ERROR(contrast_experiment_framework(
        "max_strength", CONTRAST_EXPERIMENT_FN { return max_out_strength(graph.adj_list, max_n_seeds); }));
  }
  if (params.with_pagerank) {
    RFL_RETURN_ON_ERROR(contrast_experiment_framework(
        "pagerank", CONTRAST_EXPERIMENT_FN {
          auto pagerank_params = pagerank_params_common;
          pagerank_params.transpose = false;
          return min_pagerank(graph, pagerank_params);
        }));
  }
  if (params.with_pagerank_transpose) {
    RFL_RETURN_ON_ERROR(contrast_experiment_framework(
        "pagerank_transpose", CONTRAST_EXPERIMENT_FN {
          auto pagerank_params = pagerank_params_common;
          pagerank_params.transpose = true;
          return max_pagerank(graph, pagerank_params);
        }));
  }
  if (params.with_imrank) {
    RFL_RETURN_ON_ERROR(contrast_experiment_framework(
        "imrank", CONTRAST_EXPERIMENT_FN {
          auto imrank_params = IMRankParams{
              .k = max_n_seeds,
              .n_iterations = params.imrank_n_iterations,
              .n_iterations_before_topk_fixed = *params.imrank_n_iterations_before_topk_fixed,
          };
          return imrank(graph.inv_adj_list, imrank_params);
        }));
  }
  if (params.with_rr_sketch) {
    auto& json_exp_cur = json_root["rr_sketching"];
    auto& json_time_cur = json_time_used["rr_sketching"];
    auto wim_params = WIMParams{
        .n_sketches = params.n_sketches,
        .n_seeds = params.n_seeds,
        .simulation_try_count = params.simulation_try_count,
    };
    ELOG_INFO << "==== Starts RR-sketching algorithm ====";
    RFL_RETURN_ON_ERROR(
        do_wim_experiment_get_seeds( //
            *graph_for_algorithm, *params.common, wim_params, &json_exp_cur["rr_sketch"], &json_time_cur)
            .and_then([&](VertexListList seeds_selected) -> rfl::Result<VertexListList> {
              return do_coarsening_experiment_expand_seeds( //
                  graph, coarsen_results, std::move(seeds_selected), params.n_sketches, *params.n_fast_expanding_levels,
                  *params.expanding, &json_exp_cur["expanding_seeds"], &json_time_cur["expand_seeds"]);
            })
            .transform([&](VertexListList expanded_seeds) {
              do_wim_experiment_simulate(graph, expanded_seeds, wim_params, //
                                         &json_exp_cur["simulation"], &json_time_cur["simulate"]);
              return RESULT_VOID_SUCCESS;
            }));
  }
  // Done all
  return std::move(json_root);
}

#undef CONTRAST_EXPERIMENT_FN

template <class ParamsType, class DoExperimentFn>
auto experiment_framework(int argc, char** argv, DoExperimentFn&& do_experiment_fn) -> ResultVoid try {
  auto timer = nw::util::seconds_timer{};

  return ParamsType::parse_from_args(argc, argv).and_then([&](ParamsType params) {
    init_easylog(*params.common);
    ELOGFMT(INFO, "Parameters: {:4}", params);
    auto json_fout_ptr = create_json_fout_ptr(params.common->json_output_file);

    timer.start();
    return read_directed_wim_adjacency_lists(params.common->input_file)
        .and_then([&](AdjacencyListPair<WIMEdge> read_result) {
          auto read_graph_time = timer.lap();
          auto n = graph::num_vertices(read_result.adj_list);
          auto m = read_result.adj_list.num_edges();
          ELOGFMT(INFO, "Done reading graph. |V| = {}, |E| = {}, time usage = {:.3} sec.", n, m, read_graph_time);

          // During experiment: do_experiment_fn(const AdjacencyListPair<E>&, const ParamsType&) -> rfl::Result<json>
          return std::invoke(do_experiment_fn, read_result, params).transform([&](json json_root) {
            auto json_root_str = json_root.dump(4);
            dump_to_json_fout_str(json_fout_ptr.get(), json_root_str);
            return RESULT_VOID_SUCCESS;
          });
        });
  });
}
RFL_RESULT_CATCH_HANDLER()
} // namespace

auto wim_experiment(int argc, char** argv) noexcept -> ResultVoid {
  return experiment_framework<WIMExperimentParams>(argc, argv, do_wim_experiment);
}

auto wim_coarsening_experiment(int argc, char** argv) noexcept -> ResultVoid {
  return experiment_framework<WIMCoarseningExperimentParams>(argc, argv, do_wim_coarsening_experiment);
}

auto wim_contrast_experiment(int argc, char** argv) noexcept -> ResultVoid {
  return experiment_framework<WIMContrastExperimentParams>(argc, argv, do_wim_contrast_experiment);
}
