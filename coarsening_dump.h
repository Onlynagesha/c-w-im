#pragma once

#include "coarsening_details_types.h"
#include <fmt/ranges.h>

inline auto dump_p_internal(const CoarsenedVertexDetailsBase::PInternalContainer& p_internal, //
                            size_t n_members, int indent = 0, int level = 0) -> std::string {
  indent = std::max(indent, 0);
  level = std::max(level, 0);

  auto m = n_members;
  auto rows = range(m) | TRANSFORM_VIEW(fmt::format("{::.4f}", p_internal[_1] | views::take(m)));
  if (indent <= 0) {
    return fmt::format("{}", rows);
  }
  auto res = "["s;
  for (auto i : range(m)) {
    res += '\n' + std::string(indent * (level + 1), ' ');
    res += rows[i];
  }
  res += '\n' + std::string(indent * level, ' ') + ']';
  return res;
}

inline auto dump_p_cross(const CoarsenedEdgeDetailsBase::PCrossContainer& p_cross, //
                         size_t n_members_left, size_t n_members_right,            //
                         int indent = 0, int level = 0) -> std::string {
  indent = std::max(indent, 0);
  level = std::max(level, 0);

  auto n_rows = n_members_left;
  auto n_cols = n_members_right;
  auto rows = range(n_rows) | TRANSFORM_VIEW(fmt::format("{::.4f}", p_cross[_1] | views::take(n_cols)));
  if (indent <= 0) {
    return fmt::format("{}", rows);
  }
  auto res = "["s;
  for (auto i : range(n_rows)) {
    res += '\n' + std::string(indent * (level + 1), ' ');
    res += rows[i];
  }
  res += '\n' + std::string(indent * level, ' ') + ']';
  return res;
};
