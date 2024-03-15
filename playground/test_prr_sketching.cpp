#include "dump.h"
#include "playground/sample_graph.h"
#include "utils/easylog.h"
#include "wim.h"
#include <fmt/ranges.h>

int main() {
  easylog::set_min_severity(easylog::Severity::DEBUG);

  auto [graph, inv_graph] = make_sample_wbim_graph_1();
  auto [n, m] = graph_n_m(graph);
  auto seeds = VertexSet(n, {1, 9});

  auto vertex_weights = [n]() {
    auto view = views::iota(10_vid, 10_vid + n);
    return std::vector<vertex_weight_t>(view.begin(), view.end());
  }();

  constexpr auto N_PRR_SKETCHES = 100'000;
  constexpr auto SIMULATION_TRY_COUNT = 100'000;

  auto sketches = PRRSketchSet(&graph, &inv_graph, vertex_weights, &seeds);
  rand_engine.seed(1); // Fixed seed
  sketches.append(N_PRR_SKETCHES, n);
  // ELOGFMT(INFO, "After {} attempts of PRR-sketching: sketches = {:4}", N_PRR_SKETCHES, sketches);

  easylog::set_min_severity(easylog::Severity::DEBUG);

  auto sim_res_0 = wbim_simulate_w(graph, vertex_weights, seeds, {n, {}}, SIMULATION_TRY_COUNT);

  auto selected = sketches.select(4);
  auto sim_res = wbim_simulate_w(graph, vertex_weights, seeds, {n, selected}, SIMULATION_TRY_COUNT);
  ELOGFMT(INFO, "Selected boost-vertices = {}. Simulation result = {:.4f}", selected, *sim_res - *sim_res_0);

  auto selected_by_critical = sketches.select_by_critical(4);
  auto sim_res_c = wbim_simulate_w(graph, vertex_weights, seeds, {n, selected_by_critical}, SIMULATION_TRY_COUNT);
  ELOGFMT(INFO, "Selected boost-vertices by critical = {}. Simulation result = {:.4f}", selected_by_critical,
          *sim_res_c - *sim_res_0);

  return 0;
}
