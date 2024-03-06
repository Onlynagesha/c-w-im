#include "coarsening.h"
#include "dump.h"
#include <fmt/ranges.h>
#include <nwgraph/adaptors/edge_range.hpp>

// ---- Dump functions ----

#define DEFINE_FREE_DUMP_FUNCTION_FOR_COARSENING_DETAILS(Type)                \
  auto dump(const Type& obj, int indent, int level) noexcept -> std::string { \
    return obj.dump(indent, level);                                           \
  }

COARSENING_DETAILS_TYPES(DEFINE_FREE_DUMP_FUNCTION_FOR_COARSENING_DETAILS)
#undef DEFINE_FREE_DUMP_FUNCTION_FOR_COARSENING_DETAILS

namespace {
template <ranges::forward_range Range, class ToStringFn>
  requires(std::is_invocable_r_v<std::string, ToStringFn, ranges::range_value_t<Range>>)
auto dump_array_as_braced_list(Range&& values, ToStringFn&& to_str_fn, int indent, int level) -> std::string {
  constexpr auto N_VALUES_PER_ROW = 10;
  if (indent <= 0) {
    auto joined = values | views::transform(to_str_fn) | views::join_with(", "s) | views::common;
    return '{' + std::string(joined.begin(), joined.end()) + '}';
  }
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
auto dump_array_as_braced_list(Range&& values, int width, int indent, int level) -> std::string {
  auto to_str_fn = LAMBDA_1(fmt::format("{0:{1}}", _1, width));
  return dump_array_as_braced_list(std::forward<Range>(values), to_str_fn, indent, level);
}

template <forward_range_of<std::string> ComponentRange>
auto merge_dumped_components(ComponentRange&& components, int indent, int level) -> std::string {
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

template <is_edge_property E, int IsInv>
auto dump_graph_generic(const graph::adjacency<IsInv, E>& graph, int indent, int level) -> std::string {
  auto edge_range = graph::make_edge_range<0>(remove_const(graph));
  auto components = make_reserved_vector<std::string>(graph.num_edges());
  for (auto [u, v, w] : edge_range) {
    components.push_back(fmt::format("({}, {}): {}", u, v, w));
  }
  return merge_dumped_components(components, indent, level);
}

template <class CoarseningDetailsType>
auto dump_impl(const CoarsenGraphResult<CoarseningDetailsType>& result_obj, int indent, int level) -> std::string {
  constexpr auto DECIMAL_DIGITS = 4;
  auto weights_width = // 1 : One position for the decimal point '.'
      static_cast<int>(std::log10(ranges::max(result_obj.coarsened_vertex_weights))) + DECIMAL_DIGITS + 1;
  auto weight_to_str = LAMBDA_1(fmt::format("{0:{1}.{2}f}", _1, weights_width, DECIMAL_DIGITS));
  auto weights_str = dump_array_as_braced_list(result_obj.coarsened_vertex_weights, weight_to_str, indent, level + 1);
  auto components = {
      fmt::format(".coarsened_graph = {}", dump_graph_generic(result_obj.coarsened_graph, indent, level + 1)),
      fmt::format(".coarsened_inv_graph = {}", dump_graph_generic(result_obj.coarsened_inv_graph, indent, level + 1)),
      fmt::format(".coarsened_vertex_weights = {}", weights_str),
      fmt::format(".details = {}", result_obj.details.dump(indent, level + 1))};
  return merge_dumped_components(components, indent, level);
}
} // namespace

auto CoarsenedVertexBrief::dump(int indent, int level) const noexcept -> std::string {
  indent = std::max(indent, 0);
  level = std::max(level, 0);
  auto components = {fmt::format(".members = {}", members), fmt::format(".best_seed_index = {}", best_seed_index)};
  return merge_dumped_components(components, indent, level);
}

auto CoarsenedVertexDetails::dump(int indent, int level) const noexcept -> std::string {
  indent = std::max(indent, 0);
  level = std::max(level, 0);
  auto dump_p_internal = [&](const PInternalContainer& p_internal) {
    auto m = n_members();
    auto rows = range(m) | views::transform(LAMBDA_1(fmt::format("{::.4f}", p_internal[_1] | views::take(m))));
    if (indent <= 0) {
      return fmt::format("{}", rows);
    }
    auto res = "["s;
    for (auto i : range(m)) {
      res += '\n' + std::string(indent * (level + 2), ' ');
      res += rows[i];
    }
    res += '\n' + std::string(indent * (level + 1), ' ') + ']';
    return res;
  };
  auto components = {fmt::format(".members = {}", members),
                     fmt::format(".vertex_weights = {::.4f}", vertex_weights),
                     fmt::format(".heuristics_in = {::.4f}", heuristics_in),
                     fmt::format(".heuristics_out = {::.4f}", heuristics_out),
                     fmt::format(".heuristics_out_seed = {::.4f}", heuristics_out_seed),
                     fmt::format(".p_internal = {}", dump_p_internal(p_internal)),
                     fmt::format(".p_seed_internal = {}", dump_p_internal(p_seed_internal))};
  return merge_dumped_components(components, indent, level);
}

auto CoarsenedEdgeDetails::dump(int indent, int level) const noexcept -> std::string {
  indent = std::max(indent, 0);
  level = std::max(level, 0);
  auto dump_p_cross = [&](const PCrossContainer& p_cross) {
    auto n_rows = n_members_left;
    auto n_cols = n_members_right;
    auto rows = range(n_rows) | views::transform(LAMBDA_1(fmt::format("{::.4f}", p_cross[_1] | views::take(n_cols))));
    if (indent <= 0) {
      return fmt::format("{}", rows);
    }
    auto res = "["s;
    for (auto i : range(n_rows)) {
      res += '\n' + std::string(indent * (level + 2), ' ');
      res += rows[i];
    }
    res += '\n' + std::string(indent * (level + 1), ' ') + ']';
    return res;
  };
  auto components = {fmt::format(".p_cross = {}", dump_p_cross(p_cross)),           //
                     fmt::format(".p_seed_cross = {}", dump_p_cross(p_seed_cross)), //
                     fmt::format(".merged = {}", merged)};
  return merge_dumped_components(components, indent, level);
}

auto CoarseningBrief::dump(int indent, int level) const noexcept -> std::string {
  indent = std::max(indent, 0);
  level = std::max(level, 0);
  auto groups_str = [&]() {
    auto to_group_str = LAMBDA_1(fmt::format("[{}] = {}", _1, ::dump(groups[_1])));
    return merge_dumped_components(range(n_coarsened) | views::transform(to_group_str), indent, level + 1);
  }();
  auto components = {fmt::format(".n = {}", n), fmt::format(".n_coarsened = {}", n_coarsened),
                     fmt::format(".groups = {}", groups_str)};
  return merge_dumped_components(components, indent, level);
}

auto CoarseningDetails::dump(int indent, int level) const noexcept -> std::string {
  indent = std::max(indent, 0);
  level = std::max(level, 0);
  auto groups_str = [&]() {
    auto to_group_str = LAMBDA_1(fmt::format("[{}] = {}", _1, groups[_1].dump(indent, level + 2)));
    return merge_dumped_components(range(n_coarsened) | views::transform(to_group_str), indent, level + 1);
  }();
  auto group_id_width = static_cast<int>(std::log10(n_coarsened)) + 1;
  auto components = {
      fmt::format(".n = {}", n), fmt::format(".n_coarsened = {}", n_coarsened),
      fmt::format(".group_id = {}", dump_array_as_braced_list(group_id, group_id_width, indent, level + 1)),
      fmt::format(".index_in_group = {}", dump_array_as_braced_list(index_in_group, 1, indent, level + 1)),
      fmt::format(".groups = {}", groups_str)};
  return merge_dumped_components(components, indent, level);
}

template <same_as_either<CoarseningBrief, CoarseningDetails> DetailsType>
auto CoarsenGraphResult<DetailsType>::dump(int indent, int level) const noexcept -> std::string {
  return dump_impl(*this, indent, level);
}
