#pragma once

#include "graph_types.h"
#include "utils/utils.h"
#include <fmt/format.h>
#include <nwgraph/adaptors/edge_range.hpp>

namespace dump_utils {
template <ranges::forward_range Range, class ToStringFn>
  requires(std::is_invocable_r_v<std::string, ToStringFn, ranges::range_value_t<Range>>)
inline auto dump_array_as_braced_list(Range&& values, ToStringFn&& to_str_fn, int indent, int level) -> std::string {
  if (indent <= 0) {
    auto joined = values | views::transform(to_str_fn) | views::join_with(", "s) | views::common;
    return '{' + std::string(joined.begin(), joined.end()) + '}';
  }

  constexpr auto N_VALUES_PER_ROW = 10;
  auto res = "{"s;
  auto has_value = false; // Ensures that the component range is traversed only once
  for (auto [i, x] : values | views::enumerate) {
    if (i != 0) {
      res += ", ";
    }
    if (i % N_VALUES_PER_ROW == 0) {
      res += '\n' + std::string(indent * (level + 1), ' ');
    }
    res += to_str_fn(x);
    has_value = true;
  }
  if (has_value) {
    res += '\n' + std::string(indent * level, ' ') + '}';
  } else {
    res += '}'; // Empty list: formatted as "{}"
  }
  return res;
}

template <ranges::forward_range Range>
inline auto dump_array_as_braced_list(Range&& values, int width, int indent, int level) -> std::string {
  auto to_str_fn = LAMBDA_1(fmt::format("{0:{1}}", _1, width));
  return dump_array_as_braced_list(std::forward<Range>(values), to_str_fn, indent, level);
}

// Multi-pass of the given range is required
template <ranges::forward_range Range>
  requires(std::is_integral_v<ranges::range_value_t<Range>>)
inline auto dump_integer_array_as_braced_list(Range&& values, int indent, int level) -> std::string {
  // At most 3 passes: (1) Range size checking; (2) Min-max calculation; (3) To-string operations
  if (ranges::empty(values)) {
    return "{}"; // Corner case: minmax of empty range is undefined behavior.
  }
  auto [min, max] = ranges::minmax(values);
  auto width = 1;
  for (auto x : {min, max}) {
    if (x < 0) {
      width = std::max(width, 2 + static_cast<int>(std::log10(-x))); // Plus the length of '-'
    } else if (x > 0) {
      width = std::max(width, 1 + static_cast<int>(std::log10(x)));
    }
  }
  return dump_array_as_braced_list(values, width, indent, level);
}

template <ranges::forward_range Range>
  requires(std::is_floating_point_v<ranges::range_value_t<Range>>)
inline auto dump_floating_point_array_as_braced_list(Range&& values, int indent, int level) -> std::string {
  using ValueType = ranges::range_value_t<Range>;
  if (ranges::empty(values)) {
    return "{}"; // Corner case
  }
  auto min_value = std::numeric_limits<ValueType>::max();    // MAX initially
  auto max_value = std::numeric_limits<ValueType>::lowest(); // MIN initially
  auto all_nan_or_inf = true;
  for (auto x : values) {
    if (std::isnan(x) || std::isinf(x)) {
      continue;
    }
    all_nan_or_inf = false;
    min_value = std::min(min_value, x);
    max_value = std::max(max_value, x);
  }
  if (all_nan_or_inf) {
    return dump_array_as_braced_list(values, 4, indent, level); // 4 : Max length of "nan", "inf" and "-inf"
  }

  constexpr auto MAX_DIGITS_BEFFORE_DECIMAL_POINT = 5;
  constexpr auto DECIMAL_FIXED_DIGITS = 4;

  // [0.0001, 99999.9999]
  auto threshold_min = std::pow(10.0, -DECIMAL_FIXED_DIGITS);
  auto threshold_max = std::pow(10.0, MAX_DIGITS_BEFFORE_DECIMAL_POINT) - threshold_min;

  auto width_before_decimal_point = 1;
  auto uses_exponent_representation = false;

  for (auto x : {min_value, max_value}) {
    if (x == 0.0) {
      continue;
    }
    auto has_minus_sign = (x < 0.0);
    x = std::fabs(x);
    if (x < threshold_min || x > threshold_max) {
      uses_exponent_representation = true;
    } else {
      width_before_decimal_point =
          std::max(width_before_decimal_point, has_minus_sign + 1 + static_cast<int>(std::log10(x)));
    }
  }
  if (uses_exponent_representation) {
    auto width = 6 + DECIMAL_FIXED_DIGITS;
    auto to_str_fn = LAMBDA_1(fmt::format("{0:{1}.{2}e}", _1, width, DECIMAL_FIXED_DIGITS));
    return dump_array_as_braced_list(values, to_str_fn, indent, level);
  } else {
    auto width = width_before_decimal_point + 1 + DECIMAL_FIXED_DIGITS; // 1 : Length of the decimal point '.'
    auto to_str_fn = LAMBDA_1(fmt::format("{0:{1}.{2}f}", _1, width, DECIMAL_FIXED_DIGITS));
    return dump_array_as_braced_list(values, to_str_fn, indent, level);
  }
}

template <forward_range_of<std::string> ComponentRange>
inline auto merge_dumped_components(ComponentRange&& components, int indent, int level) -> std::string {
  if (indent <= 0) {
    auto joined = components | views::join_with(", "s) | views::common;
    return '{' + std::string(joined.begin(), joined.end()) + '}';
  }
  auto res = "{"s;
  auto has_component = false; // Ensures that the component range is traversed only once
  for (auto [i, c] : components | views::enumerate) {
    if (i != 0) {
      res += ',';
    }
    res += '\n' + std::string(indent * (level + 1), ' ');
    res += c;
    has_component = true;
  }
  if (has_component) {
    res += '\n' + std::string(indent * level, ' ') + '}';
  } else {
    res += '}'; // Empty list: formatted as "{}"
  }
  return res;
}

template <class Property, int IsInv>
  requires(fmt::is_formattable<Property>::value)
inline auto dump_graph(const graph::adjacency<IsInv, Property>& graph, int indent, int level) -> std::string {
  auto edge_range = graph::make_edge_range<0>(as_non_const(graph));
  auto components = make_reserved_vector<std::string>(graph.num_edges());
  for (auto [u, v, w] : edge_range) {
    components.push_back(fmt::format("({}, {}): {}", u, v, w));
  }
  return merge_dumped_components(components, indent, level);
}
} // namespace dump_utils
