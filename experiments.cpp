#include "experiments.h"
#include "dump.h"
#include "experiments_internal_states.h"
#include "graph_connectivity.h"
#include "utils/easylog.h"
#include "utils/histogram.h"
#include "wim.h"
#include <fmt/ranges.h>
#include <magic_enum.hpp>

namespace {
using exp_states::VertexList;
using exp_states::VertexListList;

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
auto dump_vertices_selected(const AdjacencyListPair<E>& graph, std::span<const vertex_id_t> vertices_selected)
    -> std::string {
  auto [n, m] = graph.graph_n_m();
  auto res = fmt::format("Vertices selected, with |E| / |V| = {:.3f}:", 1.0 * m / n);
  for (auto seed : vertices_selected) {
    BOOST_ASSERT_MSG(seed < n, "Vertex index out of range [0, n).");
    res += fmt::format("\n\tvertex-id = {}, in-degree = {}, out-degree = {}", //
                       seed, graph.in_degree(seed), graph.out_degree(seed));
  }
  return res;
}

template <is_edge_property E>
auto write_graph_basic_information(json& json_root, const AdjacencyList<E>& graph) {
  auto n = graph.num_vertices()[0];
  auto m = graph.num_edges();
  json_root["n"] = n;
  json_root["m"] = m;
  return std::tuple{n, m};
}

template <is_edge_property E>
auto write_graph_connectivity_information(json& json_root, const AdjacencyListPair<E>& graph) {
  // # of SCCs is absent since the recursive algorithm may fail due to DFS stack limit
  auto n_wcc = n_weakly_connected_components(graph.adj_list, graph.inv_adj_list);
  json_root["n_wcc"] = n_wcc;
  return n_wcc;
}

auto do_wim_experiment_get_seeds(const WIMAdjacencyListPair& graph, const CommonExperimentParams& common,
                                 const WIMSketchingParams& sketching)
    -> rfl::Result<exp_states::WIMSketchingGetSeedsResult> {
  auto timer = nw::util::seconds_timer{};

  const auto& [adj_list, inv_adj_list, vertex_weights] = graph;
  auto [n, m] = graph_n_m(adj_list);
  auto max_n_seeds = ranges::max(sketching.n_seeds);

  auto rr_sketches_total_time_used = 0.0;
  auto rr_sketches = RRSketchSet{&inv_adj_list, vertex_weights};

  auto n_sketching_groups = sketching.n_sketches.size();
  auto res = exp_states::WIMSketchingGetSeedsResult{
      .selected_seeds = make_reserved_vector<VertexList>(n_sketching_groups),
      .sketching_info = make_reserved_vector<exp_states::WIMSketchingInfo>(n_sketching_groups),
  };
  for (auto [i, r] : sketching.n_sketches | views::enumerate) {
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
            rr_sketches.ratio_of_single_vertex_sketch() * 100.0,    // {5}
            rr_sketches.rr_sketch_total_size_str()                  // {6}
    );
    // Histogram of the distribution of RR-sketch size
    MYLOG_FMT_DEBUG("Histogram of RR-sketch sizes:\n{}",
                    make_histogram(rr_sketches.sketch_sizes(), common.histogram_shape()));

    timer.start();
    auto selected_seeds = rr_sketches.select(max_n_seeds);
    timer.stop();
    auto seed_selecting_time_usage = timer.elapsed();
    ELOGFMT(INFO, "Done selecting {} seeds with {} RR-sketches. Time used = {:.3f} sec.", //
            max_n_seeds, r, seed_selecting_time_usage);
    MYLOG_DEBUG(dump_vertices_selected(graph, selected_seeds));

    // Result
    auto info_item = exp_states::WIMSketchingInfo{
        .n_sketches = r,
        .selected_seeds = std::move(selected_seeds),
        .average_sketch_size = avg_sketch_size,
        .sketching_total_time_usage = rr_sketches_total_time_used,
        .seed_selecting_time_usage = seed_selecting_time_usage,
    };
    // A COPY into res.selected_seeds
    res.selected_seeds.push_back(info_item.selected_seeds);
    res.sketching_info.push_back(std::move(info_item));
  }
  return std::move(res);
}

auto do_wbim_experiment_get_boosted(const WBIMAdjacencyListPair& graph, const VertexSet& seeds,
                                    const CommonExperimentParams& common, const WBIMSketchingParams& sketching)
    -> rfl::Result<exp_states::WBIMSketchingGetBoostedResult> {
  auto timer = nw::util::seconds_timer{};

  const auto& [adj_list, inv_adj_list, vertex_weights] = graph;
  auto [n, m] = graph_n_m(adj_list);
  auto max_n_boosted = ranges::max(sketching.n_boosted);

  auto prr_sketches_total_time_used = 0.0;
  auto prr_sketches = PRRSketchSet{&adj_list, &inv_adj_list, vertex_weights, &seeds};

  auto n_sketching_groups = sketching.n_sketches.size();
  auto res = exp_states::WBIMSketchingGetBoostedResult{
      .selected_boosted = make_reserved_vector<VertexList>(n_sketching_groups),
      .sketching_info = make_reserved_vector<exp_states::WBIMSketchingInfo>(n_sketching_groups),
  };
  for (auto [i, r] : sketching.n_sketches | views::enumerate) {
    auto r_new = r - prr_sketches.n_sketches();
    timer.start();
    prr_sketches.append(r_new, max_n_boosted);
    timer.stop();

    auto prr_sketches_new_time_used = timer.elapsed();
    prr_sketches_total_time_used += prr_sketches_new_time_used;

    auto avg_sketch_size = prr_sketches.average_sketch_size();
    auto success_rate = prr_sketches.success_rate();
    // Log of time usage & statistics of RR-sketches
    constexpr auto msg_pattern_1 = //
        "Done appending new {1} PRR-sketches in {3:.3f} sec."
        "\n\tAverage size of all the {0} RR-sketches = {4:.3f}, success rate = {5:.3f}%)."
        "\n\tTotal time used = {2:.3f} sec.";
    ELOGFMT(INFO, msg_pattern_1,                                      //
            r, r_new,                                                 // {0}, {1}
            prr_sketches_total_time_used, prr_sketches_new_time_used, // {2}, {3}
            avg_sketch_size,                                          // {4}
            success_rate * 100.0                                      // {5}
    );
    // Histogram of the distribution of RR-sketch size
    MYLOG_FMT_DEBUG("Histogram of PRR-sketch sizes:\n{}",
                    make_histogram(prr_sketches.sketch_sizes(), common.histogram_shape()));

    timer.start();
    auto selected_boosted = prr_sketches.select(max_n_boosted);
    timer.stop();
    auto boosted_selecting_time_usage = timer.elapsed();
    ELOGFMT(INFO, "Done selecting {} boosted vertices with {} PRR-sketches. Time used = {:.3f} sec.", //
            max_n_boosted, r, boosted_selecting_time_usage);
    MYLOG_DEBUG(dump_vertices_selected(graph, selected_boosted));

    // Result
    auto info_item = exp_states::WBIMSketchingInfo{
        .n_sketches = r,
        .selected_boosted = std::move(selected_boosted),
        .average_sketch_size = avg_sketch_size,
        .sketching_success_rate = success_rate,
        .sketching_total_time_usage = prr_sketches_total_time_used,
        .boosted_selecting_time_usage = boosted_selecting_time_usage,
    };
    // A COPY into res.selected_boosted
    res.selected_boosted.push_back(info_item.selected_boosted);
    res.sketching_info.push_back(std::move(info_item));
  }
  return std::move(res);
}

auto do_wim_simulate(const WIMAdjacencyListPair& graph, std::span<const vertex_id_t> seeds_selected,
                     uint64_t simulation_try_count) -> rfl::Result<exp_states::WIMSimulationInfo> {
  const auto& [adj_list, inv_adj_list, vertex_weights] = graph;
  auto n = graph::num_vertices(adj_list);
  auto timer = nw::util::seconds_timer{};
  timer.start();
  return wim_simulate_w(adj_list, vertex_weights, {n, seeds_selected}, simulation_try_count)
      .transform([&](double sim_res) {
        timer.stop();
        auto time_used = timer.elapsed();
        ELOGFMT(INFO, "Done WIM simulation. Result = {:.3f}. Time used = {:.3f} sec.", sim_res, time_used);
        return exp_states::WIMSimulationInfo{
            .n_seeds = static_cast<vertex_id_t>(seeds_selected.size()),
            .objective_function = sim_res,
            .time_usage = time_used,
        };
      });
}

auto do_wbim_simulate(const WBIMAdjacencyListPair& graph, const VertexSet& seeds,
                      std::span<const vertex_id_t> boosted_selected, uint64_t simulation_try_count,
                      const double* base_objective_function) -> rfl::Result<exp_states::WBIMSimulationInfo> {
  const auto& [adj_list, inv_adj_list, vertex_weights] = graph;
  auto n = graph::num_vertices(adj_list);
  auto timer = nw::util::seconds_timer{};

  auto total_time_used = 0.0;
  auto base_objective_function_value = [&]() -> rfl::Result<double> {
    if (base_objective_function != nullptr) {
      return *base_objective_function;
    }
    // Base = F(EmptySet ; S), i.e. no boosted vertices involved
    timer.start();
    auto res = wbim_simulate_w(adj_list, vertex_weights, seeds, {n, {}}, simulation_try_count);
    timer.stop();
    total_time_used += timer.elapsed();
    return res;
  }();
  RFL_RETURN_ON_ERROR(base_objective_function_value);

  timer.start();
  return wbim_simulate_w(adj_list, vertex_weights, seeds, {n, boosted_selected}, simulation_try_count)
      .transform([&](double sim_res) {
        timer.stop();
        auto objective_function = sim_res - *base_objective_function_value;
        total_time_used += timer.elapsed();
        ELOGFMT(INFO, "Done WBIM simulation. Result = {:.3f} = {:.3f} = {:.3f}. Time used = {:.3f} sec.", //
                sim_res, *base_objective_function_value, objective_function, total_time_used);
        return exp_states::WBIMSimulationInfo{
            .n_boosted = static_cast<vertex_id_t>(boosted_selected.size()),
            .objective_function = objective_function,
            .time_usage = total_time_used,
        };
      });
}

auto do_wim_simulate(const WIMAdjacencyListPair& graph, std::span<const vertex_id_t> seeds_selected,
                     std::span<const vertex_id_t> n_seeds, uint64_t simulation_try_count)
    -> rfl::Result<std::vector<exp_states::WIMSimulationInfo>> {
  BOOST_ASSERT_MSG(ranges::max(n_seeds) <= ranges::size(seeds_selected), "No enough seeds!");

  auto res = make_reserved_vector<exp_states::WIMSimulationInfo>(n_seeds.size());
  // Seeds are sorted in descending order by the estimated influence
  for (auto [i, s] : n_seeds | views::enumerate) {
    auto item = do_wim_simulate(graph, seeds_selected.first(s), simulation_try_count);
    RFL_RETURN_ON_ERROR(item);
    res.push_back(std::move(*item));
  }
  return {std::move(res)};
}

auto do_wbim_simulate(const WBIMAdjacencyListPair& graph, const VertexSet& seeds,
                      std::span<const vertex_id_t> boosted_selected, std::span<const vertex_id_t> n_boosted,
                      uint64_t simulation_try_count, const double* base_objective_function)
    -> rfl::Result<std::vector<exp_states::WBIMSimulationInfo>> {
  BOOST_ASSERT_MSG(ranges::max(n_boosted) <= ranges::size(boosted_selected), "No enough boosted vertices!");

  const auto& [adj_list, inv_adj_list, vertex_weights] = graph;
  auto [n, m] = graph_n_m(adj_list);

  auto base_value = [&]() -> rfl::Result<double> {
    if (base_objective_function != nullptr) {
      return *base_objective_function;
    }
    return wbim_simulate_w(adj_list, vertex_weights, seeds, {n, {}}, simulation_try_count);
  }();
  RFL_RETURN_ON_ERROR(base_value);

  auto res = make_reserved_vector<exp_states::WBIMSimulationInfo>(n_boosted.size());
  // Boosted vertices are sorted in descending order by estimated influence
  for (auto [i, b] : n_boosted | views::enumerate) {
    auto item = do_wbim_simulate( //
        graph, seeds, boosted_selected.first(b), simulation_try_count, std::addressof(*base_value));
    RFL_RETURN_ON_ERROR(item);
    res.push_back(std::move(*item));
  }
  return {std::move(res)};
}

auto do_wim_experiment_simulate(const WIMAdjacencyListPair& graph, std::span<const VertexList> seeds_selected,
                                const WIMSketchingParams& params, uint64_t simulation_try_count)
    -> rfl::Result<std::vector<exp_states::WIMSketchingSimulationResult>> {
  BOOST_ASSERT_MSG(seeds_selected.size() == params.n_sketches.size(), "Mismatch with # of RR-sketch groups.");

  auto res = make_reserved_vector<exp_states::WIMSketchingSimulationResult>(params.n_sketches.size());
  for (auto [i, r] : params.n_sketches | views::enumerate) {
    auto items = do_wim_simulate(graph, seeds_selected[i], params.n_seeds, simulation_try_count);
    RFL_RETURN_ON_ERROR(items);
    res.push_back({.n_sketches = r, .simulation_results = std::move(*items)});
  }
  return {std::move(res)};
}

auto do_wbim_experiment_simulate(const WBIMAdjacencyListPair& graph, const VertexSet& seeds,
                                 std::span<const VertexList> boosted_selected, const WBIMSketchingParams& params,
                                 uint64_t simulation_try_count, const double* base_objective_function)
    -> rfl::Result<std::vector<exp_states::WBIMSketchingSimulationResult>> {
  BOOST_ASSERT_MSG(boosted_selected.size() == params.n_sketches.size(), "Mismatch with # of PRR-sketch groups.");

  const auto& [adj_list, inv_adj_list, vertex_weights] = graph;
  auto [n, m] = graph_n_m(adj_list);

  auto base_value = [&]() -> rfl::Result<double> {
    if (base_objective_function != nullptr) {
      return *base_objective_function;
    }
    return wbim_simulate_w(adj_list, vertex_weights, seeds, {n, {}}, simulation_try_count);
  }();
  RFL_RETURN_ON_ERROR(base_value);

  auto res = make_reserved_vector<exp_states::WBIMSketchingSimulationResult>(params.n_sketches.size());
  for (auto [i, r] : params.n_sketches | views::enumerate) {
    auto items = do_wbim_simulate( //
        graph, seeds, boosted_selected[i], params.n_boosted, simulation_try_count, std::addressof(*base_value));
    RFL_RETURN_ON_ERROR(items);
    res.push_back({.n_sketches = r, .simulation_results = std::move(*items)});
  }
  return {std::move(res)};
}

template <is_edge_property E, class CoarsenFn>
  requires(std::is_invocable_r_v<rfl::Result<typename CoarseningDataTraits<E>::BriefResult>, //
                                 CoarsenFn, const AdjacencyListPair<E>&>)
auto do_continue_coarsening_graph_impl(const AdjacencyListPair<E>& original_graph,
                                       std::vector<typename CoarseningDataTraits<E>::BriefResult>& coarsen_results,
                                       vertex_id_t destination_level, vertex_id_t coarsening_threshold,
                                       CoarsenFn&& coarsen_fn) -> rfl::Result<exp_states::CoarseningInfo> {
  using BriefResultType = typename CoarseningDataTraits<E>::BriefResult;

  auto initial_level = static_cast<vertex_id_t>(coarsen_results.size());
  if (initial_level >= destination_level) {
    // No need to coarsen the graph. Dummy value is returned.
    return exp_states::CoarseningInfo{.n_coarsened = 0, .m_coarsened = 0, .time_usage = 0.0};
  }
  auto n_coarsened = (initial_level == 0) ? original_graph.n_vertices() : coarsen_results.back().coarsened.n_vertices();

  auto timer = nw::util::seconds_timer{};
  timer.start();
  for (auto level = initial_level + 1; n_coarsened > coarsening_threshold && level <= destination_level; level++) {
    auto cur_result = [&]() -> rfl::Result<const BriefResultType*> {
      auto item = (level <= 1) ? coarsen_fn(original_graph) : coarsen_fn(coarsen_results.back().coarsened);
      return item.transform([&](const auto& item) { return &coarsen_results.emplace_back(std::move(item)); });
    }();
    RFL_RETURN_ON_ERROR(cur_result);
    n_coarsened = (*cur_result)->details.n_coarsened;
    ELOGFMT(INFO, "Finishes coarsening level {}: |V|, |E| = {}", level, (*cur_result)->coarsened.graph_n_m());
  }
  auto time_used = timer.lap();
  ELOGFMT(INFO, "Finishes coarsening from level {} to {}: Time used = {:.3f} sec.", //
          initial_level, destination_level, time_used);
  // Done coarsening to the destination level
  auto [nc, mc] = coarsen_results.back().coarsened.graph_n_m();
  return exp_states::CoarseningInfo{.n_coarsened = nc, .m_coarsened = mc, .time_usage = time_used};
}

auto do_continue_coarsening_wim_graph(const WIMAdjacencyListPair& original_graph,
                                      std::vector<WIMCoarsenGraphBriefResult>& coarsen_results,
                                      vertex_id_t destination_level, vertex_id_t coarsening_threshold,
                                      const CoarseningParams& params) -> rfl::Result<exp_states::CoarseningInfo> {
  return do_continue_coarsening_graph_impl(
      original_graph, coarsen_results, destination_level, coarsening_threshold,
      [&params](const WIMAdjacencyListPair& p) { return coarsen_wim_graph_p(p, params); });
}

auto do_continue_coarsening_wbim_graph(const WBIMAdjacencyListPair& original_graph, const VertexSet& seeds,
                                       std::vector<WBIMCoarsenGraphBriefResult>& coarsen_results,
                                       vertex_id_t destination_level, vertex_id_t coarsening_threshold,
                                       const CoarseningParams& params) -> rfl::Result<exp_states::CoarseningInfo> {
  const auto* seeds_used = coarsen_results.empty() ? &seeds : &coarsen_results.back().coarsened_seeds;
  return do_continue_coarsening_graph_impl(
      original_graph, coarsen_results, destination_level, coarsening_threshold,
      [seeds_used, &params](const WBIMAdjacencyListPair& p) { return coarsen_wbim_graph_p(p, *seeds_used, params); });
}

auto do_continue_coarsening_wim_graph(const WIMAdjacencyListPair& original_graph,
                                      std::vector<WIMCoarsenGraphBriefResult>& coarsen_results,
                                      vertex_id_t destination_level, const MultiLevelParams& multi_level_params) {
  return do_continue_coarsening_wim_graph(                //
      original_graph, coarsen_results, destination_level, //
      *multi_level_params.coarsening_threshold, *multi_level_params.coarsening);
}

auto do_continue_coarsening_wbim_graph(const WBIMAdjacencyListPair& original_graph, const VertexSet& seeds,
                                       std::vector<WBIMCoarsenGraphBriefResult>& coarsen_results,
                                       vertex_id_t destination_level, const MultiLevelParams& multi_level_params) {
  return do_continue_coarsening_wbim_graph(                      //
      original_graph, seeds, coarsen_results, destination_level, //
      *multi_level_params.coarsening_threshold, *multi_level_params.coarsening);
}

auto do_coarsen_wim_graph(const WIMAdjacencyListPair& graph, vertex_id_t destination_level,
                          vertex_id_t coarsening_threshold, const CoarseningParams& params)
    -> rfl::Result<exp_states::CoarseningResult<WIMEdge>> {
  auto res = exp_states::CoarseningResult<WIMEdge>{};
  res.coarsen_results.reserve(destination_level);
  return do_continue_coarsening_wim_graph( //
             graph, res.coarsen_results, destination_level, coarsening_threshold, params)
      .transform([&](exp_states::CoarseningInfo info) {
        res.info = info;
        return std::move(res);
      });
}

auto do_coarsen_wbim_graph(const WBIMAdjacencyListPair& graph, const VertexSet& seeds, vertex_id_t destination_level,
                           vertex_id_t coarsening_threshold, const CoarseningParams& params)
    -> rfl::Result<exp_states::CoarseningResult<WBIMEdge>> {
  auto res = exp_states::CoarseningResult<WBIMEdge>{};
  res.coarsen_results.reserve(destination_level);
  return do_continue_coarsening_wbim_graph( //
             graph, seeds, res.coarsen_results, destination_level, coarsening_threshold, params)
      .transform([&](exp_states::CoarseningInfo info) {
        res.info = info;
        return std::move(res);
      });
}

auto do_coarsen_wim_graph(const WIMAdjacencyListPair& graph, vertex_id_t destination_level,
                          const MultiLevelParams& multi_level_params) {
  return do_coarsen_wim_graph( //
      graph, destination_level, *multi_level_params.coarsening_threshold, *multi_level_params.coarsening);
}

auto do_coarsen_wbim_graph(const WBIMAdjacencyListPair& graph, const VertexSet& seeds, vertex_id_t destination_level,
                           const MultiLevelParams& multi_level_params) {
  return do_coarsen_wbim_graph( //
      graph, seeds, destination_level, *multi_level_params.coarsening_threshold, *multi_level_params.coarsening);
}

template <is_edge_property E>
auto do_expand_vertices(const AdjacencyListPair<E>& original_graph,
                        std::span<const typename CoarseningDataTraits<E>::BriefResult> coarsen_results,
                        VertexList vertices, vertex_id_t n_fast_expanding_levels, const ExpandingParams& raw_params)
    -> rfl::Result<typename exp_states::ExpansionTraits<E>::InfoType> {
  auto level = static_cast<vertex_id_t>(coarsen_results.size());
  auto rule_local = ExpandingParams{.vertex_expanding_rule = VertexExpandingRule::LOCAL};
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
              magic_enum::enum_name(raw_params.vertex_expanding_rule), back_level);
      return &raw_params;
    }();
    auto expanded_vertices = [&]() {
      if (back_level == 0) {
        if constexpr (std::is_same_v<E, WIMEdge>) {
          return expand_wim_seed_vertices(adj_list, vertex_weights, coarsen_results[0], vertices, *params_to_use)
              .transform([&](ExpandSeedResult res) { return std::move(res.expanded_seeds); });
        } else if constexpr (std::is_same_v<E, WBIMEdge>) {
          return expand_wbim_boosted_vertices(adj_list, vertex_weights, coarsen_results[0], vertices, *params_to_use)
              .transform([&](ExpandBoostedResult res) { return std::move(res.expanded_boosted); });
        } else {
          static_assert(rfl::always_false_v<E>, "Unsupported edge type.");
        }
      }
      auto& prev = coarsen_results[back_level - 1];
      auto& cur = coarsen_results[back_level];
      if constexpr (std::is_same_v<E, WIMEdge>) {
        return further_expand_wim_seed_vertices(prev, cur, vertices, *params_to_use)
            .transform([&](ExpandSeedResult res) { return std::move(res.expanded_seeds); });
      } else if constexpr (std::is_same_v<E, WBIMEdge>) {
        return further_expand_wbim_boosted_vertices(prev, cur, vertices, *params_to_use)
            .transform([&](ExpandBoostedResult res) { return std::move(res.expanded_boosted); });
      } else {
        static_assert(rfl::always_false_v<E>, "Unsupported edge type.");
      }
    }();
    RFL_RETURN_ON_ERROR(expanded_vertices);
    MYLOG_FMT_DEBUG("Seeds expansion back to level {} -> {}", back_level, expanded_vertices);
    vertices = *expanded_vertices; // Updates seeds or boosted by the expanded
  }
  timer.stop();
  auto time_used = timer.elapsed();

  if constexpr (std::is_same_v<E, WIMEdge>) {
    return exp_states::WIMExpansionInfo{.expanded_seeds = std::move(vertices), .time_usage = time_used};
  } else if constexpr (std::is_same_v<E, WBIMEdge>) {
    return exp_states::WBIMExpansionInfo{.expanded_boosted = std::move(vertices), .time_usage = time_used};
  } else {
    static_assert(rfl::always_false_v<E>, "Invalid edge type.");
  }
}

template <is_edge_property E>
auto do_expand_vertices(const AdjacencyListPair<E>& original_graph,
                        std::span<const typename CoarseningDataTraits<E>::BriefResult> coarsen_results,
                        VertexListList vertex_lists, vertex_id_t n_fast_expanding_levels,
                        const ExpandingParams& raw_params)
    -> rfl::Result<typename exp_states::ExpansionTraits<E>::ResultType> {
  using InfoType = typename exp_states::ExpansionTraits<E>::InfoType;
  using ResultType = typename exp_states::ExpansionTraits<E>::ResultType;

  auto n_groups = vertex_lists.size();
  auto res = ResultType{};
  auto [info_expanded_list, res_expanded_list_list] = [&]() {
    if constexpr (std::is_same_v<E, WIMEdge>) {
      return std::tuple{&InfoType::expanded_seeds, &ResultType::expanded_seeds};
    } else if constexpr (std::is_same_v<E, WBIMEdge>) {
      return std::tuple{&InfoType::expanded_boosted, &ResultType::expanded_boosted};
    } else {
      static_assert(rfl::always_false_v<E>, "Invalid edge type.");
    }
  }();

  (res.*res_expanded_list_list).reserve(n_groups);
  res.info_items.reserve(n_groups);

  for (auto [i, cur_seed_list] : vertex_lists | views::enumerate) {
    auto cur_expanded = do_expand_vertices( //
        original_graph, coarsen_results, cur_seed_list, n_fast_expanding_levels, raw_params);
    RFL_RETURN_ON_ERROR(cur_expanded);
    // A copy into res.expanded_vertices
    (res.*res_expanded_list_list).push_back((*cur_expanded).*info_expanded_list);
    res.info_items.push_back(std::move(*cur_expanded));
  }
  return {std::move(res)};
}

template <is_edge_property E>
auto do_expand_vertices(const AdjacencyListPair<E>& original_graph,
                        std::span<const typename CoarseningDataTraits<E>::BriefResult> coarsen_results,
                        VertexListList vertex_lists, const MultiLevelParams& multi_level_params) {
  return do_expand_vertices(original_graph, coarsen_results, std::move(vertex_lists), //
                            *multi_level_params.n_fast_expanding_levels, *multi_level_params.expanding);
}

template <class WIMParamsType>
  requires(requires(const WIMParamsType& params) {
    { *params.common } -> std::convertible_to<CommonExperimentParams>;
    { *params.sketching } -> std::convertible_to<WIMSketchingParams>;
  })
auto do_wim_experiment(const WIMAdjacencyListPair& graph, const WIMParamsType& params, json* json_root) -> ResultVoid {
  BOOST_ASSERT(json_root != nullptr);
  write_graph_basic_information(*json_root, graph.adj_list);

  return do_wim_experiment_get_seeds(graph, *params.common, *params.sketching)
      .and_then([&](exp_states::WIMSketchingGetSeedsResult get_seeds_res) {
        return do_wim_experiment_simulate( //
                   graph, get_seeds_res.selected_seeds, *params.sketching, *params.common->simulation_try_count)
            .transform([&](std::vector<exp_states::WIMSketchingSimulationResult> sim_res) {
              BOOST_ASSERT(sim_res.size() == params.sketching->n_sketches.size());
              (*json_root)["sketching"] = exp_states::to_json(get_seeds_res);
              (*json_root)["simulating"] = exp_states::to_json(sim_res);
              return RESULT_VOID_SUCCESS;
            });
      });
}

template <class WBIMParamsType>
  requires(requires(const WBIMParamsType& params) {
    { *params.common } -> std::convertible_to<CommonExperimentParams>;
    { *params.sketching } -> std::convertible_to<WBIMSketchingParams>;
  })
auto do_wbim_experiment(const WBIMAdjacencyListPair& graph, const VertexSet& seeds, const WBIMParamsType& params,
                        json* json_root) -> ResultVoid {
  BOOST_ASSERT(json_root != nullptr);
  write_graph_basic_information(*json_root, graph.adj_list);

  return do_wbim_experiment_get_boosted(graph, seeds, *params.common, *params.sketching)
      .and_then([&](exp_states::WBIMSketchingGetBoostedResult get_boosted_res) {
        auto [n, m] = graph.graph_n_m();
        auto try_count = *params.common->simulation_try_count;
        return wbim_simulate_w(graph.adj_list, graph.vertex_weights, seeds, {n, {}}, try_count)
            .and_then([&](double base_F) {
              return do_wbim_experiment_simulate( //
                  graph, seeds, get_boosted_res.selected_boosted, *params.sketching, try_count, &base_F);
            })
            .transform([&](std::vector<exp_states::WBIMSketchingSimulationResult> sim_res) {
              BOOST_ASSERT(sim_res.size() == params.sketching->n_sketches.size());
              (*json_root)["sketching"] = exp_states::to_json(get_boosted_res);
              (*json_root)["simulating"] = exp_states::to_json(sim_res);
              return RESULT_VOID_SUCCESS;
            });
      });
}

auto do_wim_sketching_3_step_process(const AdjacencyListPair<WIMEdge>& original_graph,
                                     std::span<const WIMCoarsenGraphBriefResult> coarsen_results,
                                     const CommonExperimentParams& common_params,
                                     const WIMSketchingParams& sketching_params,
                                     const MultiLevelParams& multi_level_params, json* json_root) -> ResultVoid {
  BOOST_ASSERT(json_root != nullptr);
  const auto* graph_for_algo = coarsen_results.empty() ? &original_graph : &coarsen_results.back().coarsened;
  return do_wim_experiment_get_seeds(*graph_for_algo, common_params, sketching_params)
      .and_then([&](exp_states::WIMSketchingGetSeedsResult coarsened_res) {
        (*json_root)["sketching"] = exp_states::to_json(coarsened_res);
        return do_expand_vertices(original_graph, coarsen_results, coarsened_res.selected_seeds, multi_level_params);
      })
      .and_then([&](exp_states::WIMExpansionResult expansion_res) {
        (*json_root)["expanding"] = exp_states::to_json(expansion_res);
        return do_wim_experiment_simulate( //
            original_graph, expansion_res.expanded_seeds, sketching_params, *common_params.simulation_try_count);
      })
      .transform([&](std::vector<exp_states::WIMSketchingSimulationResult> sim_res) {
        (*json_root)["simulating"] = exp_states::to_json(sim_res);
        return RESULT_VOID_SUCCESS;
      });
}

auto do_wim_sketching_3_step_process(const AdjacencyListPair<WIMEdge>& original_graph,
                                     std::span<const WIMCoarsenGraphBriefResult> coarsen_results,
                                     const WIMCoarseningExperimentParams& params, json* json_root) -> ResultVoid {
  return do_wim_sketching_3_step_process( //
      original_graph, coarsen_results, *params.common, *params.sketching, *params.multi_level, json_root);
}

auto do_wbim_sketching_3_step_process(const AdjacencyListPair<WBIMEdge>& original_graph, const VertexSet& seeds,
                                      std::span<const WBIMCoarsenGraphBriefResult> coarsen_results,
                                      const CommonExperimentParams& common_params,
                                      const WBIMSketchingParams& sketching_params,
                                      const MultiLevelParams& multi_level_params, json* json_root) -> ResultVoid {
  BOOST_ASSERT(json_root != nullptr);

  const auto* graph_for_sketching = //
      coarsen_results.empty() ? &original_graph : &coarsen_results.back().coarsened;
  const auto* seeds_for_sketching = //
      coarsen_results.empty() ? &seeds : &coarsen_results.back().coarsened_seeds;

  return do_wbim_experiment_get_boosted(*graph_for_sketching, *seeds_for_sketching, common_params, sketching_params)
      .and_then([&](exp_states::WBIMSketchingGetBoostedResult coarsened_res) {
        (*json_root)["sketching"] = exp_states::to_json(coarsened_res);
        return do_expand_vertices(original_graph, coarsen_results, coarsened_res.selected_boosted, multi_level_params);
      })
      .and_then([&](exp_states::WBIMExpansionResult expansion_res) {
        (*json_root)["expanding"] = exp_states::to_json(expansion_res);
        auto try_count = *common_params.simulation_try_count;
        return do_wbim_experiment_simulate( //
            original_graph, seeds, expansion_res.expanded_boosted, sketching_params, try_count, nullptr);
      })
      .transform([&](std::vector<exp_states::WBIMSketchingSimulationResult> sim_res) {
        (*json_root)["simulating"] = exp_states::to_json(sim_res);
        return RESULT_VOID_SUCCESS;
      });
}

auto do_wbim_sketching_3_step_process(const AdjacencyListPair<WBIMEdge>& original_graph, const VertexSet& seeds,
                                      std::span<const WBIMCoarsenGraphBriefResult> coarsen_results,
                                      const WBIMCoarseningExperimentParams& params, json* json_root) -> ResultVoid {
  return do_wbim_sketching_3_step_process( //
      original_graph, seeds, coarsen_results, *params.common, *params.sketching, *params.multi_level, json_root);
}

auto do_wim_coarsening_experiment(const AdjacencyListPair<WIMEdge>& graph, const WIMCoarseningExperimentParams& params,
                                  json* json_root) -> ResultVoid {
  BOOST_ASSERT(json_root != nullptr);
  auto timer = nw::util::seconds_timer{};
  auto& json_exp = (*json_root)["experiments"];

  const auto& [adj_list, inv_adj_list, vertex_weights] = graph;
  auto [n, m] = write_graph_basic_information(*json_root, adj_list);

  if (n <= *params.multi_level->coarsening_threshold) {
    constexpr auto msg_pattern = "Experiments fails to proceed since |V| <= threshold, with |V| = {}, threshold = {}";
    return rfl::Error{fmt::format(msg_pattern, n, *params.multi_level->coarsening_threshold)};
  }
  // Step 1: Level 0, i.e. initial graph
  {
    ELOG_INFO << "Starts Level 0: solves WIM problem on the original graph.";
    auto& json_exp_cur = json_exp.emplace_back(json{{"level", 0}});
    auto n_wcc = write_graph_connectivity_information(json_exp_cur, graph);
    ELOGFMT(INFO, "# of weakly connected components = {}", n_wcc);
    RFL_RETURN_ON_ERROR(do_wim_experiment(graph, params, &json_exp_cur));
  }

  // Step 2: Coarsening until the size threshold
  auto coarsen_results = std::vector<WIMCoarsenGraphBriefResult>{};
  for (auto level = 1_vid, n_coarsened = n; n_coarsened > params.multi_level->coarsening_threshold; ++level) {
    auto& json_exp_cur = json_exp.emplace_back(json{{"level", level}});
    // Step 2.1: Gets the coarsened graph
    auto coarsen_success = //
        do_continue_coarsening_wim_graph(graph, coarsen_results, level, *params.multi_level)
            .transform([&](exp_states::CoarseningInfo info) {
              json_exp_cur["coarsening"] = exp_states::to_json(info);
              return RESULT_VOID_SUCCESS;
            });
    RFL_RETURN_ON_ERROR(coarsen_success);
    const auto& cur_coarsen_result = coarsen_results.back();
    n_coarsened = cur_coarsen_result.coarsened.n_vertices();

    // Checks connectivity for debugging & data analysis
    auto n_wcc = write_graph_connectivity_information(json_exp_cur, cur_coarsen_result.coarsened);
    ELOGFMT(INFO, "# of weakly connected components = {}", n_wcc);
    // Step 2.2: Solves on the coarsened graph -> exp_states::anding seeds -> Simulation to estimate F(S)
    RFL_RETURN_ON_ERROR(do_wim_sketching_3_step_process(graph, coarsen_results, params, &json_exp_cur));
  }
  return RESULT_VOID_SUCCESS;
}

auto do_wbim_coarsening_experiment(const AdjacencyListPair<WBIMEdge>& graph, const VertexSet& seeds,
                                   const WBIMCoarseningExperimentParams& params, json* json_root) -> ResultVoid {
  BOOST_ASSERT(json_root != nullptr);
  auto timer = nw::util::seconds_timer{};
  auto& json_exp = (*json_root)["experiments"];

  const auto& [adj_list, inv_adj_list, vertex_weights] = graph;
  auto [n, m] = write_graph_basic_information(*json_root, adj_list);

  if (n <= *params.multi_level->coarsening_threshold) {
    constexpr auto msg_pattern = "Experiments fails to proceed since |V| <= threshold, with |V| = {}, threshold = {}";
    return rfl::Error{fmt::format(msg_pattern, n, *params.multi_level->coarsening_threshold)};
  }
  // Step 1: Level 0, i.e. initial graph
  {
    ELOG_INFO << "Starts Level 0: solves WBIM problem on the original graph.";
    auto& json_exp_cur = json_exp.emplace_back(json{{"level", 0}});
    auto n_wcc = write_graph_connectivity_information(json_exp_cur, graph);
    ELOGFMT(INFO, "# of weakly connected components = {}", n_wcc);
    RFL_RETURN_ON_ERROR(do_wbim_experiment(graph, seeds, params, &json_exp_cur));
  }

  // Step 2: Coarsening until the size threshold
  auto coarsen_results = std::vector<WBIMCoarsenGraphBriefResult>{};
  for (auto level = 1_vid, n_coarsened = n; n_coarsened > params.multi_level->coarsening_threshold; ++level) {
    auto& json_exp_cur = json_exp.emplace_back(json{{"level", level}});
    // Step 2.1: Gets the coarsened graph
    auto coarsen_success = //
        do_continue_coarsening_wbim_graph(graph, seeds, coarsen_results, level, *params.multi_level)
            .transform([&](exp_states::CoarseningInfo info) {
              json_exp_cur["coarsening"] = exp_states::to_json(info);
              return RESULT_VOID_SUCCESS;
            });
    RFL_RETURN_ON_ERROR(coarsen_success);
    const auto& cur_coarsen_result = coarsen_results.back();
    n_coarsened = cur_coarsen_result.coarsened.n_vertices();

    // Checks connectivity for debugging & data analysis
    auto n_wcc = write_graph_connectivity_information(json_exp_cur, cur_coarsen_result.coarsened);
    ELOGFMT(INFO, "# of weakly connected components = {}", n_wcc);
    // Step 2.2: Solves on the coarsened graph -> exp_states::anding seeds -> Simulation to estimate F(S)
    RFL_RETURN_ON_ERROR(do_wbim_sketching_3_step_process(graph, seeds, coarsen_results, params, &json_exp_cur));
  }
  return RESULT_VOID_SUCCESS;
}

auto do_wim_contrast_experiment(const AdjacencyListPair<WIMEdge>& graph, const WIMContrastExperimentParams& params,
                                json* json_root) -> ResultVoid {
  auto timer = nw::util::seconds_timer{};
  auto [n, m] = write_graph_basic_information(*json_root, graph.adj_list);

  auto coarsen_results = //
      do_coarsen_wim_graph(graph, params.coarsening_level, *params.multi_level)
          .transform([&](exp_states::CoarseningResult<WIMEdge> results) {
            (*json_root)["coarsening"] = exp_states::to_json(results.info);
            return std::move(results.coarsen_results);
          });
  RFL_RETURN_ON_ERROR(coarsen_results);

  if (params.with_rr_sketch) {
    auto& json_exp_cur = (*json_root)["rr_sketching"];
    auto sketching_params = WIMSketchingParams{
        .n_sketches = params.n_sketches,
        .n_seeds = params.n_seeds,
    };
    ELOG_INFO << "==== Starts RR-sketching algorithm ====";
    RFL_RETURN_ON_ERROR(do_wim_sketching_3_step_process( //
        graph, *coarsen_results, *params.common, sketching_params, *params.multi_level, &json_exp_cur));
  }

#define CONTRAST_EXPERIMENT_FN [&](const WIMAdjacencyListPair& graph) -> rfl::Result<VertexList>

  auto contrast_experiment_framework = [&](std::string_view algorithm_name, auto&& algorithm_fn) -> ResultVoid {
    ELOGFMT(INFO, "==== Starts contrast algorithm: {} ====", algorithm_name);
    auto& json_exp_cur = (*json_root)[algorithm_name];
    timer.start();
    // Step 1: Performs the algorithm
    const auto* graph_for_algo = coarsen_results->empty() ? &graph : &coarsen_results->back().coarsened;
    return std::invoke(algorithm_fn, *graph_for_algo)
        .and_then([&](VertexList seeds) {
          timer.stop();
          auto algo_time_used = timer.elapsed();
          json_exp_cur["algorithm"] = {
              {"seeds_before_expanding", seeds},
              {"time_used", algo_time_used},
          };
          ELOGFMT(INFO, "Finished algorithm '{}' in {:.3f} sec.", algorithm_name, algo_time_used);
          MYLOG_FMT_DEBUG("Seed vertices before expansion = {}", seeds);
          // Step 2: exp_states::ands seed vertices
          return do_expand_vertices(graph, *coarsen_results, {std::move(seeds)}, *params.multi_level);
        })
        .and_then([&](exp_states::WIMExpansionResult expansion_result) {
          BOOST_ASSERT(expansion_result.info_items.size() == 1);
          json_exp_cur["expanding"] = exp_states::to_json(expansion_result.info_items[0]);
          // Step 3: Simulation with expanded seeds
          auto try_count = *params.common->simulation_try_count;
          return do_wim_simulate(graph, expansion_result.expanded_seeds[0], params.n_seeds, try_count);
        })
        .transform([&](std::vector<exp_states::WIMSimulationInfo> sim_info) {
          json_exp_cur["simulating"] = exp_states::to_json(sim_info);
          return RESULT_VOID_SUCCESS;
        });
  };

  auto max_n_seeds = ranges::max(params.n_seeds);
  if (params.with_max_degree) {
    RFL_RETURN_ON_ERROR(contrast_experiment_framework(
        "max_degree", CONTRAST_EXPERIMENT_FN { return max_out_degree(graph.adj_list, max_n_seeds); }));
  }
  if (params.with_max_strength) {
    RFL_RETURN_ON_ERROR(contrast_experiment_framework(
        "max_strength", CONTRAST_EXPERIMENT_FN { return max_out_strength(graph.adj_list, max_n_seeds); }));
  }
  if (params.with_pagerank) {
    auto pagerank_params = PagerankParams{
        .damping_factor = *params.pagerank_damping_factor,
        .epsilon = *params.pagerank_epsilon,
        .k = max_n_seeds,
        .n_iterations = params.pagerank_n_iterations,
        .transpose = true,
        .uses_vertex_weight = true,
        .uses_edge_weight = true,
    };
    RFL_RETURN_ON_ERROR(contrast_experiment_framework(
        "pagerank", CONTRAST_EXPERIMENT_FN { return max_pagerank(graph, pagerank_params); }));
  }
  if (params.with_imrank) {
    auto imrank_params = IMRankParams{
        .k = max_n_seeds,
        .n_iterations = params.imrank_n_iterations,
        .n_iterations_before_topk_fixed = *params.imrank_n_iterations_before_topk_fixed,
    };
    RFL_RETURN_ON_ERROR(contrast_experiment_framework(
        "imrank", CONTRAST_EXPERIMENT_FN { return imrank(graph.inv_adj_list, imrank_params); }));
  }
  // Done all
  return RESULT_VOID_SUCCESS;
}

#undef CONTRAST_EXPERIMENT_FN

template <class ParamsType, class DoExperimentFn>
auto wim_experiment_framework(int argc, char** argv, DoExperimentFn&& do_experiment_fn) -> ResultVoid try {
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

          // During experiment: do_experiment_fn(const AdjacencyListPair<E>&, const ParamsType&) -> ResultVoid
          auto json_root = json{};
          return std::invoke(do_experiment_fn, read_result, params, &json_root).transform([&](rfl::Nothing) {
            auto json_root_str = json_root.dump(4);
            dump_to_json_fout_str(json_fout_ptr.get(), json_root_str);
            return RESULT_VOID_SUCCESS;
          });
        });
  });
}
RFL_RESULT_CATCH_HANDLER()

template <class ParamsType, class DoExperimentFn>
auto wbim_experiment_framework(int argc, char** argv, DoExperimentFn&& do_experiment_fn) -> ResultVoid try {
  auto timer = nw::util::seconds_timer{};

  return ParamsType::parse_from_args(argc, argv).and_then([&](ParamsType params) {
    init_easylog(*params.common);
    ELOGFMT(INFO, "Parameters: {:4}", params);
    auto json_fout_ptr = create_json_fout_ptr(params.common->json_output_file);

    timer.start();
    return read_directed_wbim_adjacency_lists(params.common->input_file)
        .and_then([&](AdjacencyListPair<WBIMEdge> read_result) -> ResultVoid {
          auto [n, m] = read_result.graph_n_m();
          if (*params.n_seeds_to_generate <= 0 || *params.n_seeds_to_generate >= n) {
            constexpr auto msg_pattern = "Invalid # of seeds to generate (input: {}): "
                                         "Must be in the range [1, |V|-1] (with |V| = {})";
            return rfl::Error{fmt::format(msg_pattern, *params.n_seeds_to_generate, n - 1)};
          }
          auto seeds = VertexSet{n, max_out_strength(read_result.adj_list, *params.n_seeds_to_generate)};
          ELOGFMT(INFO, "Seeds generated for WBIM experiment:\n{}", //
                  dump_vertices_selected(read_result, seeds.vertex_list));

          timer.stop();
          auto read_graph_time = timer.elapsed();
          ELOGFMT(INFO, "Done reading graph. |V| = {}, |E| = {}, time usage = {:.3} sec.", n, m, read_graph_time);

          // During experiment: do_experiment_fn(const AdjacencyListPair<E>&, const ParamsType&) -> ResultVoid
          auto json_root = json{};
          return std::invoke(do_experiment_fn, read_result, seeds, params, &json_root).transform([&](rfl::Nothing) {
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
  return wim_experiment_framework<WIMSketchingExperimentParams>( //
      argc, argv, do_wim_experiment<WIMSketchingExperimentParams>);
}

auto wbim_experiment(int argc, char** argv) noexcept -> ResultVoid {
  return wbim_experiment_framework<WBIMSketchingExperimentParams>( //
      argc, argv, do_wbim_experiment<WBIMSketchingExperimentParams>);
}

auto wim_coarsening_experiment(int argc, char** argv) noexcept -> ResultVoid {
  return wim_experiment_framework<WIMCoarseningExperimentParams>(argc, argv, do_wim_coarsening_experiment);
}

auto wbim_coarsening_experiment(int argc, char** argv) noexcept -> ResultVoid {
  return wbim_experiment_framework<WBIMCoarseningExperimentParams>(argc, argv, do_wbim_coarsening_experiment);
}

auto wim_contrast_experiment(int argc, char** argv) noexcept -> ResultVoid {
  return wim_experiment_framework<WIMContrastExperimentParams>(argc, argv, do_wim_contrast_experiment);
}
