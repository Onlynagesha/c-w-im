#pragma once

#include "graph_types.h"
#include "utils/utils.h"
#include <fmt/format.h>
#include <nwgraph/adaptors/edge_range.hpp>

namespace dump_utils {
template <ranges::forward_range Range, class ToStringFn>
  requires(std::is_invocable_r_v<std::string, ToStringFn, ranges::range_value_t<Range>>)
inline auto dump_array_as_braced_list(Range&& values, ToStringFn&& to_str_fn, int indent, int level) -> std::string {
  if (ranges::empty(values)) {
    return "{}"; // Corner case for empty range
  }
  if (indent <= 0) {
    auto joined = values | views::transform(to_str_fn) | views::join_with(", "s) | views::common;
    return '{' + std::string(joined.begin(), joined.end()) + '}';
  }

  constexpr auto N_VALUES_PER_ROW = 10;
  auto res = "{"s;
  for (auto [i, x] : values | views::enumerate) {
    if (i != 0) {
      res += ", ";
    }
    if (i % N_VALUES_PER_ROW == 0) {
      res += '\n' + std::string(indent * (level + 1), ' ');
    }
    res += to_str_fn(x);
  }
  res += '\n' + std::string(indent * level, ' ') + '}';
  return res;
}

template <ranges::forward_range Range>
inline auto dump_array_as_braced_list(Range&& values, int width, int indent, int level) -> std::string {
  auto to_str_fn = LAMBDA_1(fmt::format("{0:{1}}", _1, width));
  return dump_array_as_braced_list(std::forward<Range>(values), to_str_fn, indent, level);
}

template <ranges::forward_range Range>
  requires(std::is_integral_v<ranges::range_value_t<Range>>)
inline auto dump_integer_array_as_braced_list(Range&& values, int indent, int level) -> std::string {
  if (ranges::empty(values)) {
    return "{}"; // Corner case for empty range
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
  return dump_array_as_braced_list(std::forward<Range>(values), width, indent, level);
}

template <forward_range_of<std::string> ComponentRange>
inline auto merge_dumped_components(ComponentRange&& components, int indent, int level) -> std::string {
  if (indent <= 0) {
    auto joined = components | views::join_with(", "s) | views::common;
    return '{' + std::string(joined.begin(), joined.end()) + '}';
  }
  auto res = "{"s;
  for (const auto& [i, c] : components | views::enumerate) {
    if (i != 0) {
      res += ',';
    }
    res += '\n' + std::string(indent * (level + 1), ' ');
    res += c;
  }
  res += '\n' + std::string(indent * level, ' ') + '}';
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
