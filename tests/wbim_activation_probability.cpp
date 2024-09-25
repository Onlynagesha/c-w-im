#define BOOST_TEST_MODULE "Activation Probability in WBIM Algorithm"

#include "tests/sample_graph.h"
#include "wim.h"
#include <boost/test/unit_test.hpp>
#include <fmt/ranges.h>

template <std::floating_point T>
auto check_floating_point_range(std::string_view lhs_name, std::span<const T> lhs, //
                                std::string_view rhs_name, std::span<const T> rhs, T eps = 1e-6) -> std::string {
  auto err_msg = ""s;
  if (lhs.size() != rhs.size()) {
    err_msg = fmt::format("Size mismatch: {} vs. {}", lhs.size(), rhs.size());
    goto check_floating_point_range_failure;
  }
  for (auto [i, pair] : views::zip(lhs, rhs) | views::enumerate) {
    auto [x, y] = pair;
    if (std::fabs(x - y) > eps) {
      err_msg = fmt::format("Value mismatch at index {}: {} vs. {}", i, x, y);
      goto check_floating_point_range_failure;
    }
  }
  return err_msg; // Returns empty error message as success
check_floating_point_range_failure:
  err_msg += fmt::format("\n{}\n{}\n", dump_array(lhs_name, lhs), dump_array(rhs_name, rhs));
  return err_msg; // Fails
}

template <ranges::contiguous_range LHS, ranges::contiguous_range RHS>
  requires(std::is_same_v<ranges::range_value_t<LHS>, ranges::range_value_t<RHS>> &&
           std::is_floating_point_v<ranges::range_value_t<LHS>>)
auto check_floating_point_range(std::string_view lhs_name, LHS&& lhs, //
                                std::string_view rhs_name, RHS&& rhs, //
                                ranges::range_value_t<LHS> eps = 1e-6) -> std::string {
  using T = ranges::range_value_t<LHS>;
  return check_floating_point_range(lhs_name, std::span<const T>{lhs}, rhs_name, std::span<const T>{rhs}, eps);
}

#define CHECK_FLOATING_POINT_RANGE(lhs, rhs) check_floating_point_range(#lhs, lhs, #rhs, rhs)

BOOST_AUTO_TEST_CASE(graphC) {
  auto [graph, inv_graph] = make_sample_wbim_graph("C");
  auto n = graph::num_vertices(graph);
  auto seeds = VertexSet(n, {1, 4, 6, 7, 8, 10, 12, 14, 16, 17});

  auto out = wbim_activation_probability_from_seeds(inv_graph, seeds, n);
  BOOST_CHECK(out || !"wbim_activation_probability_from_seeds() returns error.");
  if (out) {
    auto [seed_prob, seed_prob_boosted] = *out;
    // (1) seed_prob
    const auto ans = std::vector<edge_probability_t>{
        0.34140295, 0, 0.495205, 0.3691156, 0, 0.3, 0, 0, 0, 0.42679894, 0, 0.37, 0, 0.52898085, 0, 0.40999997, 0, 0,
    };
    BOOST_CHECK_EQUAL(CHECK_FLOATING_POINT_RANGE(seed_prob, ans), "");
    // (2) seed_prob_boosted
    const auto ans_boosted = std::vector<edge_probability_t>{
        0.6019707, 0, 0.80588, 0.6689191, 0, 0.6, 0, 0, 0, 0.74864817, 0, 0.74, 0, 0.8521062, 0, 0.82, 0, 0,
    };
    BOOST_CHECK_EQUAL(CHECK_FLOATING_POINT_RANGE(seed_prob_boosted, ans_boosted), "");
  }
}
