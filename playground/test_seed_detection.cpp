#include "playground/sample_graph.h"
#include "utils/easylog.h"
#include "wim.h"
#include <fmt/ranges.h>

int main() {
  easylog::set_min_severity(easylog::Severity::TRACE);

  auto [graph, inv_graph] = make_sample_wbim_graph_3();
  auto n = graph::num_vertices(graph);
  auto seeds = VertexSet(n, {1, 4, 6, 7, 8, 10, 12, 14, 16, 17});

  auto [seed_prob, seed_prob_boosted] = *wbim_detect_probability_from_seeds(inv_graph, seeds, n);
  ELOG_INFO << dump_array("seed_prob", seed_prob);
  ELOG_INFO << dump_array("seed_prob_boosted", seed_prob_boosted);

  easylog::set_min_severity(easylog::Severity::DEBUG);

  auto sim_res = simulate_p(graph, seeds, nullptr, 1'000'000);
  ELOG_INFO << dump_array("sim_res", *sim_res);

  for (auto v : range(n)) {
    if (seeds.contains(v)) {
      continue;
    }
    auto sim_res_boosted = simulate_p(graph, seeds, VertexSet(n, {v}), 1'000'000);
    ELOG_INFO << dump_array(fmt::format("sim_res_boosted[{}]", v), *sim_res_boosted);
  }

  return 0;
}

/*
Elements of array 'seed_prob' (size = 18):
        seed_prob[0] = 0.34140295
        seed_prob[2] = 0.495205
        seed_prob[3] = 0.3691156
        seed_prob[5] = 0.3
        seed_prob[9] = 0.42679894
        seed_prob[11] = 0.37
        seed_prob[13] = 0.52898085
        seed_prob[15] = 0.40999997

Elements of array 'seed_prob_boosted' (size = 18):
        seed_prob_boosted[0] = 0.6019707
        seed_prob_boosted[2] = 0.80588
        seed_prob_boosted[3] = 0.6689191
        seed_prob_boosted[5] = 0.6
        seed_prob_boosted[9] = 0.74864817
        seed_prob_boosted[11] = 0.74
        seed_prob_boosted[13] = 0.8521062
        seed_prob_boosted[15] = 0.82
*/