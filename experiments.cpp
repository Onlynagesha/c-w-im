#include "experiments.h"
#include "dump.h"
#include "utils/graph_connectivity.h"
#include "utils/histogram.h"
#include "wim.h"
#include <fmt/ranges.h>
#include <nlohmann/json.hpp>

namespace {
using SelectedSeedLists = std::vector<std::vector<vertex_id_t>>;

inline auto init_easylog(const CommonExperimentParams& params) {
  easylog::set_min_severity(params.log_severity);
  if (!params.log_output_file.empty()) {
    easylog::init_log(params.log_severity, params.log_output_file, false, params.log_console);
  }
}

// Returns nullptr if out_path is empty
inline auto create_json_fout_ptr(const std::string& out_path) -> std::unique_ptr<std::ofstream> {
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

inline auto dump_to_json_fout_str(std::ofstream* fout, std::string_view contents) -> void {
  if (fout != nullptr) {
    ELOGFMT(DEBUG, "JSON output: {}", contents);
    (*fout) << contents;
  } else {
    ELOGFMT(INFO, "JSON output: {}", contents);
  }
}

// json_root is expected to be an array, each of whose element contains a key "r" that represents the # of RR-sketches
// Example: json_root = [
//   { "r": 10000, "something": {} },
//   { "r": 50000, "something": {}, key: {} } // <-- i = 1
// ]
auto add_key_to_ith_array_element_by_r(json& json_root, std::string_view key, size_t i, size_t r) -> json& {
  BOOST_ASSERT_MSG(json_root.empty() || json_root.is_array(), "Non-array JSON object is not accepted.");
  if (i < json_root.size()) {
    // Accesses existing json_root[i] and creates a key
    json_root[i]["r"] = r;
    return json_root[i][key];
  } else {
    // Appends { "r": r, key: ... } to json_root
    while (i > json_root.size() - 1) {
      json_root.push_back({}); // Fills with empty object until i == json_root.size() - 1
    }
    json_root.push_back({{"r", r}});
    return json_root.back()[key];
  }
}

/*
json_root will be added the following contents: {
  "n": 1500, // # of vertices in the graph
  "m": 37500, // # of directed edges in the graph
}
*/
auto write_graph_basic_information(json& json_root, const AdjacencyList<WIMEdge>& graph) {
  auto n = graph.num_vertices()[0];
  auto m = graph.num_edges();
  json_root["n"] = n;
  json_root["m"] = m;
  return std::tuple{n, m};
}

/*
json_root will be added the following contents: {
  "n": 1500, // # of vertices in the graph
  "m": 37500, // # of directed edges in the graph
  "num_sketches": [10000, 100000], // # of RR-sketches to be generated
  "num_seeds": [5, 10]) // # of seed vertices to be selected
}
*/
auto write_wim_basic_information(json& json_root, const AdjacencyList<WIMEdge>& graph, const WIMParams& params) {
  auto nm_tuple = write_graph_basic_information(json_root, graph);
  json_root["num_sketches"] = params.num_sketches;
  json_root["num_seeds"] = params.num_seeds;
  return nm_tuple;
}

auto write_graph_connectivity_information(json& json_root, const AdjacencyList<WIMEdge>& graph,
                                          const InvAdjacencyList<WIMEdge>& inv_graph) {
  auto n_scc = n_strongly_connected_components(graph);
  auto n_wcc = n_weakly_connected_components(graph, inv_graph);
  ELOGFMT(INFO, "# of strongly connected components = {}; # of weakly connected components = {}", n_scc, n_wcc);
  json_root["n_scc"] = n_scc;
  json_root["n_wcc"] = n_wcc;
  return std::tuple{n_scc, n_wcc};
}

/*
Each object in the list json_experiments will be added keys "average_sketch_size" and "selected_seeds"
like the following (If json_experiments is an empty list, then new objects with keys
"r", "average_sketch_size" and "selected_seeds" will be created): [
  {
    "r": 10000,
    "average_sketch_size": 16.2, // Each RR-sketch contains 16.2 vertices in average
      // All the 10 seeds selected, sorted by estimated gain in descending order
      // e.g. The first 5 seeds, [0, 1275, 365, 45, 902], is the result when num_seeds = 5
    "selected_seeds": [ 0, 1275, 365, 45, 902, 125, 28, 692, 1234, 345 ]
  },
  {
    "r": 100000,
    "average_sketch_size": 16.5,
    "selected_seeds": [ 0, 902, 365, 1275, 45, 125, 1234, 28, 777, 345 ]
  }
]
The object json_time_used will be added the following contents: {
  "rr_sketch": [
    { "r": 10000, "seconds": 2.5 },  // Time usage of 10'000 RR-sketches in seconds
    { "r": 100000, "seconds": 24.0 } // ... 100'000 RR-sketches ...
  ],
  "select_seeds": [
    { "r": 10000, "seconds": 0.07 }, // Time usage of selecting 10 seeds with 10'000 RR-sketches
    { "r": 100000, "seconds": 0.56 } // ... with 100'000 RR-sketches
  ],
}
*/
template <class ParamsType>
auto do_wim_experiment_get_seeds(const AdjacencyList<WIMEdge>& adj_list, const InvAdjacencyList<WIMEdge>& inv_adj_list,
                                 std::span<const vertex_weight_t> vertex_weights, const ParamsType& params,
                                 json& json_experiments, json& json_time_used) noexcept
    -> rfl::Result<SelectedSeedLists> try {
  auto timer = nw::util::seconds_timer{};
  auto [n, m] = graph_n_m(adj_list);

  auto rr_sketches_total_time_used = 0.0;
  auto rr_sketches = RRSketchSet{&inv_adj_list, vertex_weights};

  auto& json_time_used_rr_sketch = json_time_used["rr_sketch"];
  auto& json_time_used_select_seeds = json_time_used["select_seeds"];

  auto res = SelectedSeedLists{};
  for (auto [i, r] : params.wim->num_sketches | views::enumerate) {
    auto r_new = r - rr_sketches.num_sketches();
    timer.start();
    rr_sketches.append(r_new);
    timer.stop();

    auto rr_sketches_new_time_used = timer.elapsed();
    rr_sketches_total_time_used += rr_sketches_new_time_used;

    auto avg_sketch_size = rr_sketches.average_sketch_size();
    // Log of time usage & statistics of RR-sketches
    constexpr auto msg_pattern_1 = //
        "Done appending new {1} RR-sketches in {3:.3f} sec. "
        "Average size of all the {0} RR-sketches = {4:.3f} "
        "(with {5:.2f}% of RR-sketches that contains 1 vertex only.); "
        "Total time used = {2:.3f} sec. ";
    ELOGFMT(INFO, msg_pattern_1, r, r_new, rr_sketches_total_time_used, rr_sketches_new_time_used, avg_sketch_size,
            rr_sketches.percentage_of_single_vertex_sketch());
    // Histogram of the distribution of RR-sketch size
    ELOGFMT(DEBUG, "Histogram of RR-sketch sizes:\n{}",
            make_histogram(rr_sketches.sketch_sizes(), params.common->histogram_shape()));
    // JSON
    add_key_to_ith_array_element_by_r(json_time_used_rr_sketch, "seconds", i, r) = rr_sketches_total_time_used;
    add_key_to_ith_array_element_by_r(json_experiments, "average_sketch_size", i, r) = avg_sketch_size;

    auto max_n_seeds = params.wim->num_seeds.back();
    timer.start();
    auto seed_list = rr_sketches.select(max_n_seeds);
    timer.stop();
    constexpr auto msg_pattern_2 = "Done selecting {} seeds with {} RR-sketches. Time used = {:.3f} sec.";
    ELOGFMT(INFO, msg_pattern_2, max_n_seeds, r, timer.elapsed());
    ELOG_DEBUG << [&] {
      auto res = fmt::format("Seeds selected (with average degree = {:.3f}):", 1.0 * m / n);
      for (auto s : seed_list) {
        auto deg = graph::degree(inv_adj_list, s) + graph::degree(adj_list, s); // In-deg + Out-deg
        res += fmt::format("\n\tid = {}, degree = {}", s, deg);
      }
      return res;
    }();
    add_key_to_ith_array_element_by_r(json_time_used_select_seeds, "seconds", i, r) = timer.elapsed();
    add_key_to_ith_array_element_by_r(json_experiments, "seeds_selected", i, r) = seed_list;
    res.push_back(std::move(seed_list));
  }
  return std::move(res);
}
RFL_RESULT_CATCH_HANDLER()

template <class ParamsType>
auto do_wim_experiment_get_seeds(const AdjacencyListPair<WIMEdge>& graph, const ParamsType& params,
                                 json& json_experiments, json& json_time_used) noexcept {
  return do_wim_experiment_get_seeds( //
      graph.adj_list, graph.inv_adj_list, graph.vertex_weights, params, json_experiments, json_time_used);
}

template <class DetailsType, class ParamsType>
auto do_wim_experiment_get_seeds(const CoarsenGraphResult<DetailsType>& graph, const ParamsType& params,
                                 json& json_experiments, json& json_time_used) noexcept {
  return do_wim_experiment_get_seeds(graph.coarsened_graph, graph.coarsened_inv_graph, graph.coarsened_vertex_weights,
                                     params, json_experiments, json_time_used);
}

/*
Each object in the list json_experiments will be added a key "simulation_result" like the following
(If json_experiments is an empty list, then new objects with keys "r" and "simulation_result" will be created): [
  { "r": 10000, "simulation_result": [
      { "s": 5, "objective_function": 632.4 },
      { "s": 10, "objective_function": 1172.5 } ]
  },
  { "r": 100000, "simulation_result": [
      { "s": 5, "objective_function": 798.4 },
      { "s": 10, "objective_function": 1375.5 } ]
  }
]
The object json_time_used will be added the following contents: {
  "simulate": [
    { "r": 10000, "by_n_seeds": [
        { "s": 5, "seconds": 1.6 },   // Time usage of simulating with the first 5 seeds with 10'000 RR-sketches
        { "s": 10, "seconds": 2.2 } ] // ... all the 10 seeds ... 10'000 RR-sketches
    },
    { "r": 100000, "by_n_seeds": [
        { "s": 5, "seconds": 1.9 },   // ... first 5 seeds ... 100'000 RR-sketches
        { "s": 10, "seconds": 2.7 } ] // ... all the 10 seeds ... 100'000 RR-sketches
    }
  ] // Adds a key "simulate" whose value is a list to json_time_used
}
*/
template <class ParamsType>
auto do_wim_experiment_simulate(const AdjacencyListPair<WIMEdge>& graph, const SelectedSeedLists& seeds_selected,
                                const ParamsType& params, json& json_experiments, json& json_time_used) noexcept
    -> ResultVoid try {
  auto timer = nw::util::seconds_timer{};
  const auto& [adj_list, inv_adj_list, vertex_weights] = graph;
  auto n = graph::num_vertices(adj_list);

  auto& json_time_used_simulate = json_time_used["simulate"];
  for (auto [i, r] : params.wim->num_sketches | views::enumerate) {
    auto& json_experiments_dest = add_key_to_ith_array_element_by_r(json_experiments, "simulation_result", i, r);
    auto& json_time_used_dest = add_key_to_ith_array_element_by_r(json_time_used_simulate, "by_n_seeds", i, r);

    for (auto s : params.wim->num_seeds) {
      auto seed_set = VertexSet{n, views::counted(seeds_selected[i].begin(), s)};
      timer.start();
      wim_simulate_w(adj_list, vertex_weights, seed_set, *params.wim->simulation_try_count)
          .and_then([&](double sim_res) -> ResultVoid {
            timer.stop();
            constexpr auto msg_pattern_2 =
                "Done simulation with the first {} seeds. Result = {:.3f}. Time used = {:.3f} sec.";
            ELOGFMT(INFO, msg_pattern_2, s, sim_res, timer.elapsed());
            json_experiments_dest.push_back({{"s", s}, {"objective_function", sim_res}});
            json_time_used_dest.push_back({{"s", s}, {"seconds", timer.elapsed()}});
            return RESULT_VOID_SUCCESS;
          })
          .or_else([&](const rfl::Error& e) -> ResultVoid {
            // Reports the error to log output and continues the experiment
            ELOGFMT(ERROR, "\nFailed simulation with {} seeds: `{}'.", s, e.what());
            return RESULT_VOID_SUCCESS;
          });
    }
  }
  return RESULT_VOID_SUCCESS;
}
RFL_RESULT_CATCH_HANDLER()

/*
Example of JSON output: {
  "n": 1500, "m": 37500, "num_sketches": [10000, 100000], "num_seeds": [5, 10]),
  "experiments": [
    { "r": 10000, "average_sketch_size": 16.2,
      "selected_seeds": [ 0, 1275, 365, 45, 902, 125, 28, 692, 1234, 345 ],
      "simulation_result": [
        { "s": 5, "objective_function": 632.4 },
        { "s": 10, "objective_function": 1172.5 },
      ]
    },
    { "r": 100000, ... }
  ],
  "time_used": {
    "rr_sketch": [ { "r": 10000, "seconds": 2.5 }, { "r": 100000, "seconds": 24.0 } ],
    "select_seeds": [ { "r": 10000, "seconds": 0.07 }, { "r": 100000, "seconds": 0.56 } ],
    "simulate": [
      { "r": 10000, "by_n_seeds": [ { "s": 5, "seconds": 1.6 }, { "s": 10, "seconds": 2.2 } ] },
      { "r": 100000, "by_n_seeds": [ { "s": 5, "seconds": 1.9 }, { "s": 10, "seconds": 2.7 } ] }
    ]
  }
}
*/
auto do_wim_experiment(const AdjacencyListPair<WIMEdge>& graph, const WIMExperimentParams& params) noexcept
    -> rfl::Result<json> try {
  auto json_root = json{};
  auto& json_experiments = json_root["experiments"];
  auto& json_time_used = json_root["time_used"];
  write_wim_basic_information(json_root, graph.adj_list, *params.wim);

  return do_wim_experiment_get_seeds(graph, params, json_experiments, json_time_used)
      .and_then([&](SelectedSeedLists seeds_selected) {
        return do_wim_experiment_simulate(graph, seeds_selected, params, json_experiments, json_time_used);
      })
      .transform([&](auto) { return std::move(json_root); });
}
RFL_RESULT_CATCH_HANDLER()

/*
Example of JSON output: {
  "n": 1500, "m": 37500, "num_sketches": [10000, 100000], "num_seeds": [5, 10]),
  "experiments": [
    { "level": 0, "n": 750, "m": 21600, "rr_sketch": [
      { "r": 10000, "average_sketch_size": 16.2, "selected_seeds": [ 0, 1275, 365, 45, 902, 125, 28, 692, 1234, 345 ] },
      { "r": 100000, "average_sketch_size": ..., "selected_seeds": [ ... ] }
    ], "simulation_result": [
      { "s": 5, "objective_function": 625.4 },
      { "s": 10, "objective_function": 1143.5 },
    ] },
    { "level": 1, "n": 375, "m": 14400, "rr_sketch": [
      { "r": 10000, "average_sketch_size": ..., "selected_seeds": [ ... ] },
      { "r": 100000, "average_sketch_size": ..., "selected_seeds": [ ... ] }
    ], "expanding_seeds": [
      { "r": 10000, "expanded_seeds": [ 0, 1275, 45, 365, 902, 28, 125, 1234, 345, 789 ] },
      { "r": 100000, "expanded_seeds": [ ... ] }
    ], "simulation_result": [
      { "s": 5, "objective_function": ... },
      { "s": 10, "objective_function": ... },
    ] },
    { "level": 2, ... }
  ],
  "time_used": [
    { "level": 0, "coarsen": 0.23,
      "rr_sketch": [ { "r": 10000, "seconds": 2.5 }, { "r": 100000, "seconds": 24.0 } ],
      "select_seeds": [ { "r": 10000, "seconds": 0.07 }, { "r": 100000, "seconds": 0.56 } ],
      "simulate": [
        { "r": 10000, "by_n_seeds": [ { "s": 5, "seconds": 1.6 }, { "s": 10, "seconds": 2.2 } ] },
        { "r": 100000, "by_n_seeds": [ { "s": 5, "seconds": 1.9 }, { "s": 10, "seconds": 2.7 } ] }
    },
    { "level": 1,
      "rr_sketch": [ { "r": 10000, "seconds": ... }, { "r": 100000, "seconds": ... } ],
      "select_seeds": [ { "r": 10000, "seconds": ... }, { "r": 100000, "seconds": ... } ],
      "expand_seeds": [ { "r": 10000, "seconds": 0.87 }, { "r": 100000, "seconds": 1.02 } ],
      "simulate": [
        { "r": 10000, "by_n_seeds": [ { "s": 5, "seconds": ... }, { "s": 10, "seconds": ... } ] },
        { "r": 100000, "by_n_seeds": [ { "s": 5, "seconds": ... }, { "s": 10, "seconds": ... } ] }
    },
  ]
}
*/
auto do_wim_coarsening_experiment(const AdjacencyListPair<WIMEdge>& graph,
                                  const WIMCoarseningExperimentParams& params) noexcept -> rfl::Result<json> try {
  auto timer = nw::util::seconds_timer{};

  auto json_root = json{};
  auto& json_experiments_root = json_root["experiments"]; // As a list
  auto& json_time_used_root = json_root["time_used"];     // As a list

  auto make_json_cur_items = [&](vertex_id_t level) {
    json_experiments_root.push_back({{"level", level}});
    json_time_used_root.push_back({{"level", level}});
    return std::tuple{&json_experiments_root.back(), &json_time_used_root.back()};
  };

  const auto& [adj_list, inv_adj_list, vertex_weights] = graph;
  auto [n, m] = write_wim_basic_information(json_root, adj_list, *params.wim);
  if (n <= params.coarsening_threshold) {
    constexpr auto msg_pattern = "Experiments fails to proceed since |V| <= threshold, with |V| = {}, threshold = {}";
    return rfl::Error{fmt::format(msg_pattern, n, *params.coarsening_threshold)};
  }
  // Step 1: Level 0, i.e. initial graph
  ELOG_INFO << "Starts Level 0: solves WIM problem on the original graph.";
  {
    auto [json_experiments, json_time_used] = make_json_cur_items(0);
    // JSON object of each level (including level 0) contains "n" and "m" that represent the corresponding graph size
    write_graph_basic_information(*json_experiments, adj_list);
    auto& json_rr_sketch = (*json_experiments)["rr_sketch"];
    // "rr_sketch": [ { "r": 10000, ... }, { "r": 50000, ... } ], see above
    auto tmp = do_wim_experiment_get_seeds(graph, params, json_rr_sketch, *json_time_used)
                   .and_then([&](SelectedSeedLists seeds) {
                     return do_wim_experiment_simulate(graph, seeds, params, json_rr_sketch, *json_time_used);
                   });
    RFL_RETURN_ON_ERROR(tmp);
  }
  // Step 2: Coarsening until the size threshold
  auto coarsen_results = std::vector<CoarsenGraphBriefResult>{};
  for (auto level = 1_vid, n_coarsened = n; n_coarsened > params.coarsening_threshold; ++level) {
    auto [json_experiments, json_time_used] = make_json_cur_items(level);
    // Step 2.1: Gets the coarsened graph
    timer.start();
    if (level == 1) {
      ELOG_INFO << "Starts graph coarsening level 1: coarsens the original graph.";
      coarsen_results.push_back(coarsen_wim_graph_w(adj_list, inv_adj_list, vertex_weights, *params.coarsening));
    } else {
      ELOGFMT(INFO, "Starts graph coarsening level {}: coarsens the previous level.", level);
      coarsen_results.push_back(further_coarsen_wim_graph(coarsen_results.back(), *params.coarsening));
    }
    timer.stop();
    const auto& cur_coarsen_result = coarsen_results.back();
    auto [n_cur, m_cur] = write_graph_basic_information(*json_experiments, cur_coarsen_result.coarsened_graph);
    ELOGFMT(INFO, "Done graph coarsening level {}: |V| => {}, |E| => {}; Time used = {:.3f} sec.", //
            level, n_cur, m_cur, timer.elapsed());
    n_coarsened = n_cur;
    (*json_time_used)["coarsen"] = timer.elapsed();
    // Step 2: Perform seed selection with the coarsened graph
    auto seed_lists =
        do_wim_experiment_get_seeds(cur_coarsen_result, params, (*json_experiments)["rr_sketch"], *json_time_used);
    RFL_RETURN_ON_ERROR(seed_lists);
    // Seed list that corresponds to each # of RR-sketches
    for (auto [i, sl] : *seed_lists | views::enumerate) {
      auto n_sketches = params.wim->num_sketches[i];
      ELOGFMT(DEBUG, "Seed list #{}: Initial seeds with r = {}: {}", i, n_sketches, sl);
      timer.start();
      // Rewinds to level 0 (i,e. seeds of the original graph)
      for (auto back_level : range(level) | views::reverse) {
        auto expanding_params = [&]() {
          if (back_level < params.n_fast_expanding_levels) {
            ELOGFMT(DEBUG, "Uses fast S_LOCAL configuration during expansion back to level {}", back_level);
            return ExpandingParams{.seed_expanding_rule = SeedExpandingRule::LOCAL};
          } else {
            ELOGFMT(DEBUG, "Uses user-provided configuration during expansion back to level {}", back_level);
            return *params.expanding;
          }
        }();
        auto expand_res = [&]() {
          if (back_level == 0) {
            return expand_wim_seed_vertices_w(adj_list, vertex_weights, coarsen_results[0], sl, expanding_params);
          } else {
            auto& last = coarsen_results[back_level - 1];
            auto& cur = coarsen_results[back_level];
            return further_expand_wim_seed_vertices(last, cur, sl, expanding_params);
          }
        }();
        RFL_RETURN_ON_ERROR(expand_res);
        ELOGFMT(DEBUG, "Seeds expansion back to level {} -> {}", back_level, expand_res->expanded_seeds);
        sl = expand_res->expanded_seeds;
      }
      timer.stop();
      ELOGFMT(DEBUG, "Finishes expansion of seed list #{} with r = {}. Time used = {:.3f} sec.", //
              i, n_sketches, timer.elapsed());
      (*json_experiments)["expanding_seeds"].push_back({{"r", n_sketches}, {"expanded_seeds", sl}});
      (*json_time_used)["expand_seeds"].push_back({{"r", n_sketches}, {"seconds", timer.elapsed()}});
    }
    // Step 3: Simulation with expanded seeds
    auto sim_res = do_wim_experiment_simulate( //
        graph, *seed_lists, params, (*json_experiments)["simulation"], *json_time_used);
    RFL_RETURN_ON_ERROR(sim_res);
  }
  return std::move(json_root);
}
RFL_RESULT_CATCH_HANDLER()

template <class ParamsType, class DoExperimentFn>
auto experiment_framework(int argc, char** argv, DoExperimentFn&& do_experiment_fn) noexcept -> ResultVoid try {
  auto timer = nw::util::seconds_timer{};

  return ParamsType::parse_from_args(argc, argv).and_then([&](ParamsType params) {
    init_easylog(*params.common);
    ELOGFMT(INFO, "Parameters: {:4}", params);
    auto json_fout_ptr = create_json_fout_ptr(params.common->json_output_file);

    timer.start();
    return read_wim_graph_data(params.common->input_file).and_then([&](AdjacencyListPair<WIMEdge> read_result) {
      auto read_graph_time = timer.lap();
      auto n = graph::num_vertices(read_result.adj_list);
      auto m = read_result.adj_list.num_edges();
      ELOGFMT(INFO, "Done reading graph. |V| = {}, |E| = {}, time usage = {:.3} sec.", n, m, read_graph_time);

      // During experiment: do_experiment_fn(const AdjacencyListPair<E>&, const ParamsType&)
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
