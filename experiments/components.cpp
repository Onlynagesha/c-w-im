#include "experiments/components.h"
#include "experiments/io.h"
#include "utils/histogram.h"
#include "wim.h"
#include <fmt/ranges.h>
#include <magic_enum.hpp>

using exp_states::VertexList;
using exp_states::VertexListList;

namespace exp_components {
// ---- Group A: Seed Generation ----

auto generate_wbim_seeds(const WBIMAdjacencyListPair& graph, const WBIMSeedGeneratingParams& params)
    -> rfl::Result<VertexSet> {
  auto [n, m] = graph.graph_n_m();
  if (**params.n_seeds >= n) {
    constexpr auto msg_pattern = "Failed to generate seeds for WBIM experiment since the required # of seeds "
                                 "(which is {}) is not less than |V| (which is {}).";
    return rfl::Error{fmt::format(msg_pattern, **params.n_seeds, n)};
  }

  auto n_candidates = static_cast<vertex_id_t>(n * (**params.candidate_ratio));
  if (n_candidates <= **params.n_seeds) {
    constexpr auto msg_pattern = "Failed to generate seeds for WBIM experiment since the # of candidates "
                                 "({} * {} = {}) is no more than the # of seeds (which {})";
    return rfl::Error{fmt::format(msg_pattern, n, **params.candidate_ratio, n_candidates, **params.n_seeds)};
  }
  MYLOG_FMT_DEBUG("Seed generation: {} candidates from {} vertices.", n_candidates, n);
  auto candidates = max_out_strength(graph.adj_list, n_candidates);

  auto n_wcc = n_weakly_connected_components(graph.adj_list, graph.inv_adj_list);
  if (n_wcc > *params.max_n_wcc_without_seeds) {
    constexpr auto msg_pattern = "Failed to generate seeds for WBIM experiment since the WCC constraint "
                                 "(<= {} WCCs required) is below the graph itself (# of WCCs = {})";
    return rfl::Error{fmt::format(msg_pattern, *params.max_n_wcc_without_seeds, n_wcc)};
  }

  auto selected_list = std::vector<vertex_id_t>(**params.n_seeds);
  auto local_rand_engine = std::mt19937{1u}; // Fixed seed
  for (auto try_index = uint64_t{0}; try_index < **params.max_try_count; try_index++) {
    // Uses local rand engine (whose seed is fixed) to ensure fixed generation result.
    ranges::sample(candidates, selected_list.begin(), selected_list.size(), local_rand_engine);
    auto n_wccs_without_selected = n_weakly_connected_components(graph.adj_list, graph.inv_adj_list, selected_list);
    if (n_wccs_without_selected > *params.max_n_wcc_without_seeds) {
      continue; // Failed this time.
    }
    constexpr auto msg_pattern = "Done seed selection after {} trials. # of WCCs after removing seeds = {}. "
                                 "Details of selected seeds:\n{}";
    ELOGFMT(INFO, msg_pattern, 1 + try_index, n_wccs_without_selected,
            exp_io::dump_vertices_selected(graph, selected_list));
    return VertexSet{n, std::move(selected_list)};
  }
  return rfl::Error{fmt::format( //
      "Failed to generate seeds for WBIM experiment after {} trials.", **params.max_try_count)};
}

// ---- Group B: Seed / Boosted Vertices Selection ----

auto do_wim_experiment_get_seeds(const WIMAdjacencyListPair& graph, const WIMSketchingParams& sketching_params,
                                 const GetSeedsParams& other_params)
    -> rfl::Result<exp_states::WIMSketchingGetSeedsResult> {
  auto timer = nw::util::seconds_timer{};

  const auto& [adj_list, inv_adj_list, vertex_weights] = graph;
  auto [n, m] = graph_n_m(adj_list);
  auto max_n_seeds = ranges::max(sketching_params.n_seeds);

  auto rr_sketches_total_time_used = 0.0;
  auto rr_sketches = RRSketchSet{&inv_adj_list, vertex_weights};

  auto n_sketching_groups = sketching_params.n_sketches.size();
  auto res = exp_states::WIMSketchingGetSeedsResult{
      .selected_seeds = make_reserved_vector<VertexList>(n_sketching_groups),
      .sketching_info = make_reserved_vector<exp_states::WIMSketchingInfo>(n_sketching_groups),
  };
  for (auto [i, r] : sketching_params.n_sketches | views::enumerate) {
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
                    make_histogram(rr_sketches.sketch_sizes(), other_params.histogram_shape));

    timer.start();
    auto selected_seeds = rr_sketches.select(max_n_seeds);
    timer.stop();
    auto seed_selecting_time_usage = timer.elapsed();
    ELOGFMT(INFO, "Done selecting {} seeds with {} RR-sketches. Time used = {:.3f} sec.", //
            max_n_seeds, r, seed_selecting_time_usage);
    MYLOG_DEBUG(exp_io::dump_vertices_selected(graph, selected_seeds));

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
                                    const WBIMSketchingParams& sketching_params, const GetBoostedParams& other_params)
    -> rfl::Result<exp_states::WBIMSketchingGetBoostedResult> {
  auto timer = nw::util::seconds_timer{};

  const auto& [adj_list, inv_adj_list, vertex_weights] = graph;
  auto [n, m] = graph_n_m(adj_list);
  auto max_n_boosted = ranges::max(sketching_params.n_boosted);

  auto prr_sketches_total_time_used = 0.0;
  auto prr_sketches = PRRSketchSet{&adj_list, &inv_adj_list, vertex_weights, &seeds};

  auto n_sketching_groups = sketching_params.n_sketches.size();
  auto res = exp_states::WBIMSketchingGetBoostedResult{
      .selected_boosted = make_reserved_vector<VertexList>(n_sketching_groups),
      .sketching_info = make_reserved_vector<exp_states::WBIMSketchingInfo>(n_sketching_groups),
  };
  for (auto [i, r] : sketching_params.n_sketches | views::enumerate) {
    auto r_new = r - prr_sketches.n_sketches();
    timer.start();
    prr_sketches.append(r_new, max_n_boosted);
    timer.stop();

    auto prr_sketches_new_time_used = timer.elapsed();
    prr_sketches_total_time_used += prr_sketches_new_time_used;

    auto avg_sketch_n_vertices = prr_sketches.average_sketch_n_vertices();
    auto avg_sketch_n_edges = prr_sketches.average_sketch_n_edges();
    auto total_sketch_size_bytes = prr_sketches.sketch_total_size_bytes();
    auto success_rate = prr_sketches.success_rate();

    // Log of time usage & statistics of RR-sketches
    constexpr auto msg_pattern_1 = //
        "Done appending new {1} PRR-sketches in {3:.3f} sec."
        "\n\tAverage |V|, |E| in PRR-sketches = {4:.3f}, {5:.3f}."
        "\n\tTotal size of PRR-sketches = {6}. success rate = {7:.3f}%."
        "\n\tTotal time used = {2:.3f} sec.";
    ELOGFMT(INFO, msg_pattern_1,                                      //
            r, r_new,                                                 // {0}, {1}
            prr_sketches_total_time_used, prr_sketches_new_time_used, // {2}, {3}
            avg_sketch_n_vertices, avg_sketch_n_edges,                // {4}, {5}
            size_bytes_to_memory_str(total_sketch_size_bytes),        // {6}
            success_rate * 100.0                                      // {7}
    );
    // Histogram of the distribution of RR-sketch size
    MYLOG_FMT_DEBUG("Distribution of |V| in PRR_sketches:\n{}",
                    make_histogram(prr_sketches.sketch_n_vertices(), other_params.histogram_shape));

    timer.start();
    auto selected_boosted = prr_sketches.select(max_n_boosted);
    timer.stop();
    auto boosted_selecting_time_usage = timer.elapsed();
    ELOGFMT(INFO, "Done selecting {} boosted vertices with {} PRR-sketches. Time used = {:.3f} sec.", //
            max_n_boosted, r, boosted_selecting_time_usage);
    MYLOG_DEBUG(exp_io::dump_vertices_selected(graph, selected_boosted));

    timer.start();
    auto selected_boosted_by_critical = prr_sketches.select_by_critical(max_n_boosted);
    timer.stop();
    auto boosted_selecting_by_critical_time_usage = timer.elapsed();
    ELOGFMT(INFO, "Done selecting {} boosted with critical points in {} PRR-sketches. Time used = {:.3f} sec.", //
            max_n_boosted, r, boosted_selecting_by_critical_time_usage);
    MYLOG_DEBUG(exp_io::dump_vertices_selected(graph, selected_boosted_by_critical));

    // Result
    auto info_item = exp_states::WBIMSketchingInfo{
        .n_sketches = r,
        .selected_boosted = std::move(selected_boosted),
        .selected_boosted_by_critical = std::move(selected_boosted_by_critical),
        .average_sketch_n_vertices = avg_sketch_n_vertices,
        .average_sketch_n_edges = avg_sketch_n_edges,
        .total_sketch_size_bytes = total_sketch_size_bytes,
        .sketching_success_rate = success_rate,
        .sketching_total_time_usage = prr_sketches_total_time_used,
        .boosted_selecting_time_usage = boosted_selecting_time_usage,
        .boosted_selecting_time_usage_by_critical = boosted_selecting_by_critical_time_usage,
    };
    // A COPY into res.selected_boosted
    res.selected_boosted.push_back(info_item.selected_boosted);
    // A COPY into res.selected_boosted_by_critical
    res.selected_boosted_by_critical.push_back(info_item.selected_boosted_by_critical);
    res.sketching_info.push_back(std::move(info_item));
  }
  return std::move(res);
}

// ---- Group C: Simulation ----

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

// ---- Group D: Graph Coarsening ----

namespace {
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
} // namespace

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

// ---- Group E: Expanding ----

namespace {
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
      // back_level > 0
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
  res.time_usage.reserve(n_groups);

  for (auto [i, cur_seed_list] : vertex_lists | views::enumerate) {
    auto cur_expanded = do_expand_vertices( //
        original_graph, coarsen_results, cur_seed_list, n_fast_expanding_levels, raw_params);
    RFL_RETURN_ON_ERROR(cur_expanded);
    // info item -> result
    (res.*res_expanded_list_list).push_back(std::move((*cur_expanded).*info_expanded_list));
    res.time_usage.push_back((*cur_expanded).time_usage);
  }
  return {std::move(res)};
}
} // namespace

auto do_expand_wim_vertices(const AdjacencyListPair<WIMEdge>& original_graph,
                            std::span<const WIMCoarsenGraphBriefResult> coarsen_results, VertexList vertices,
                            vertex_id_t n_fast_expanding_levels, const ExpandingParams& params)
    -> rfl::Result<exp_states::WIMExpansionInfo> {
  return do_expand_vertices(original_graph, coarsen_results, std::move(vertices), n_fast_expanding_levels, params);
}

auto do_expand_wim_vertex_lists(const AdjacencyListPair<WIMEdge>& original_graph,
                                std::span<const WIMCoarsenGraphBriefResult> coarsen_results,
                                VertexListList vertex_lists, vertex_id_t n_fast_expanding_levels,
                                const ExpandingParams& params) -> rfl::Result<exp_states::WIMExpansionResult> {
  return do_expand_vertices(original_graph, coarsen_results, std::move(vertex_lists), n_fast_expanding_levels, params);
}

auto do_expand_wbim_vertices(const AdjacencyListPair<WIMEdge>& original_graph,
                             std::span<const WIMCoarsenGraphBriefResult> coarsen_results, VertexList vertices,
                             vertex_id_t n_fast_expanding_levels, const ExpandingParams& params)
    -> rfl::Result<exp_states::WIMExpansionInfo> {
  return do_expand_vertices(original_graph, coarsen_results, std::move(vertices), n_fast_expanding_levels, params);
}

auto do_expand_wbim_vertex_lists(const AdjacencyListPair<WBIMEdge>& original_graph,
                                 std::span<const WBIMCoarsenGraphBriefResult> coarsen_results,
                                 VertexListList vertex_lists, vertex_id_t n_fast_expanding_levels,
                                 const ExpandingParams& params) -> rfl::Result<exp_states::WBIMExpansionResult> {
  return do_expand_vertices(original_graph, coarsen_results, std::move(vertex_lists), n_fast_expanding_levels, params);
}

} // namespace exp_components
