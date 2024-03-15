#include "coarsening.h"
#include "dump.h"
#include "dump_utils.h"
#include <fmt/ranges.h>
#include <nwgraph/adaptors/edge_range.hpp>

// ---- Dump functions ----

#define DEFINE_FREE_DUMP_FUNCTION_FOR_COARSENING_DETAILS(Type)                \
  auto dump(const Type& obj, int indent, int level) noexcept -> std::string { \
    return obj.dump(indent, level);                                           \
  }

COARSENING_DETAILS_TYPES(DEFINE_FREE_DUMP_FUNCTION_FOR_COARSENING_DETAILS)
#undef DEFINE_FREE_DUMP_FUNCTION_FOR_COARSENING_DETAILS

auto WIMCoarsenedVertexBrief::dump(int indent, int level) const noexcept -> std::string {
  indent = std::max(indent, 0);
  level = std::max(level, 0);

  auto components = {
      fmt::format(".members = {}", members),
      fmt::format(".best_seed_index = {}", best_seed_index),
  };
  return dump_utils::merge_dumped_components(components, indent, level);
}

auto WIMCoarsenedVertexDetails::dump(int indent, int level) const noexcept -> std::string {
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

  auto components = {
      fmt::format(".members = {}", members),
      fmt::format(".vertex_weights = {::.4f}", vertex_weights),
      fmt::format(".heuristics_in = {::.4f}", heuristics_in),
      fmt::format(".heuristics_out = {::.4f}", heuristics_out),
      fmt::format(".heuristics_out_seed = {::.4f}", heuristics_out_seed),
      fmt::format(".p_internal = {}", dump_p_internal(p_internal)),
      fmt::format(".p_seed_internal = {}", dump_p_internal(p_seed_internal)),
  };
  return dump_utils::merge_dumped_components(components, indent, level);
}

auto WIMCoarsenedEdgeDetails::dump(int indent, int level) const noexcept -> std::string {
  indent = std::max(indent, 0);
  level = std::max(level, 0);

  auto dump_p_cross = [&](const PCrossContainer& p_cross) {
    auto n_rows = n_members_left;
    auto n_cols = n_members_right;
    auto rows = range(n_rows) | TRANSFORM_VIEW(fmt::format("{::.4f}", p_cross[_1] | views::take(n_cols)));
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

  auto components = {
      fmt::format(".p_cross = {}", dump_p_cross(p_cross)),
      fmt::format(".p_seed_cross = {}", dump_p_cross(p_seed_cross)),
      fmt::format(".merged = {}", merged),
  };
  return dump_utils::merge_dumped_components(components, indent, level);
}

auto WIMCoarseningBrief::dump(int indent, int level) const noexcept -> std::string {
  indent = std::max(indent, 0);
  level = std::max(level, 0);

  auto groups_str = [&]() {
    auto to_group_str = LAMBDA_1(fmt::format("[{}] = {}", _1, ::dump(groups[_1])));
    return dump_utils::merge_dumped_components( //
        range(n_coarsened) | views::transform(to_group_str), indent, level + 1);
  }();

  auto components = {
      fmt::format(".n = {}", n),
      fmt::format(".n_coarsened = {}", n_coarsened),
      fmt::format(".groups = {}", groups_str),
  };
  return dump_utils::merge_dumped_components(components, indent, level);
}

auto WIMCoarseningDetails::dump(int indent, int level) const noexcept -> std::string {
  indent = std::max(indent, 0);
  level = std::max(level, 0);

  auto groups_str = [&]() {
    auto to_group_str = LAMBDA_1(fmt::format("[{}] = {}", _1, groups[_1].dump(indent, level + 2)));
    return dump_utils::merge_dumped_components( //
        range(n_coarsened) | views::transform(to_group_str), indent, level + 1);
  }();
  auto group_id_width = static_cast<int>(std::log10(n_coarsened)) + 1;

  auto components = {
      fmt::format(".n = {}", n),
      fmt::format(".n_coarsened = {}", n_coarsened),
      fmt::format(".group_id = {}", //
                  dump_utils::dump_array_as_braced_list(group_id, group_id_width, indent, level + 1)),
      fmt::format(".index_in_group = {}", //
                  dump_utils::dump_array_as_braced_list(index_in_group, 1, indent, level + 1)),
      fmt::format(".groups = {}", groups_str),
  };
  return dump_utils::merge_dumped_components(components, indent, level);
}

template <same_as_either<WIMCoarseningBrief, WIMCoarseningDetails> DetailsType>
auto WIMCoarsenGraphResult<DetailsType>::dump(int indent, int level) const noexcept -> std::string {
  constexpr auto DECIMAL_DIGITS = 4;

  auto weights_width = // 1 : One position for the decimal point '.'
      static_cast<int>(std::log10(ranges::max(coarsened.vertex_weights))) + DECIMAL_DIGITS + 1;
  auto weight_to_str = LAMBDA_1(fmt::format("{0:{1}.{2}f}", _1, weights_width, DECIMAL_DIGITS));
  auto weights_str = dump_utils::dump_array_as_braced_list( //
      coarsened.vertex_weights, weight_to_str, indent, level + 1);

  auto components = {
      fmt::format(".coarsened.adj_list = {}", //
                  dump_utils::dump_graph(coarsened.adj_list, indent, level + 1)),
      fmt::format(".coarsened.inv_adj_list = {}", //
                  dump_utils::dump_graph(coarsened.inv_adj_list, indent, level + 1)),
      fmt::format(".coarsened.vertex_weights = {}", weights_str),
      fmt::format(".details = {}", details.dump(indent, level + 1)),
  };
  return dump_utils::merge_dumped_components(components, indent, level);
}
