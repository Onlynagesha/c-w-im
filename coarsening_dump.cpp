#include "coarsening.h"
#include "dump.h"
#include "dump_utils.h"
#include <fmt/ranges.h>
#include <forward_list>
#include <nwgraph/adaptors/edge_range.hpp>

namespace {
using PInternalContainer = CoarsenedVertexDetailsBase::PInternalContainer;
using PCrossContainer = CoarsenedEdgeDetailsBase::PCrossContainer;

auto dump_p_internal(const PInternalContainer& p_internal, size_t n_members, int indent, int level) -> std::string {
  BOOST_ASSERT(level >= 0);
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

auto dump_p_cross(const PCrossContainer& p_cross, size_t n_members_left, size_t n_members_right, int indent, int level)
    -> std::string {
  BOOST_ASSERT(level >= 0);
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

auto coarsened_vertex_details_base_components(const CoarsenedVertexDetailsBase& base, int indent, int level)
    -> std::forward_list<std::string> {
  BOOST_ASSERT(level >= 0);
  return std::forward_list{
      fmt::format(".members = {}", base.members),
      fmt::format(".vertex_weights = {::.4f}", base.vertex_weights),
      fmt::format(".heuristics_in = {::.4f}", base.heuristics_in),
      fmt::format(".heuristics_out = {::.4f}", base.heuristics_out),
      fmt::format(".p_internal = {}", dump_p_internal(base.p_internal, base.n_members(), indent, level + 1)),
  };
}

auto coarsened_edge_details_base_components(const CoarsenedEdgeDetailsBase& base, int indent, int level)
    -> std::forward_list<std::string> {
  BOOST_ASSERT(level >= 0);
  auto p_cross_str = dump_p_cross(base.p_cross, base.n_members_left, base.n_members_right, indent, level + 1);
  return std::forward_list{
      fmt::format(".n_members_left = {}", base.n_members_left),
      fmt::format(".n_members_right = {}", base.n_members_right),
      fmt::format(".p_cross = {}", std::move(p_cross_str)),
  };
}

template <is_coarsened_vertex_info VertexInfoType>
auto coarsening_info_base_components(const CoarseningInfoBase<VertexInfoType>& base, int indent, int level)
    -> std::forward_list<std::string> {
  BOOST_ASSERT(level >= 0);
  auto groups_str = [&]() {
    if (indent <= 0) {
      return fmt::format("{}", base.groups);
    }
    auto components = make_reserved_vector<std::string>(base.n_coarsened);
    for (auto [i, g] : base.groups | views::enumerate) {
      components.push_back(fmt::format("[{}] = {}", i, g.dump(indent, level + 2)));
    }
    return dump_utils::merge_dumped_components(components, indent, level + 1);
  }();
  return std::forward_list{
      fmt::format(".n = {}", base.n),
      fmt::format(".n_coarsened = {}", base.n_coarsened),
      fmt::format(".groups = {}", std::move(groups_str)),
  };
}

template <is_edge_property E, class BriefOrDetailsType>
auto coarsen_graph_result_base_components(const CoarsenGraphResultBase<E, BriefOrDetailsType>& base, //
                                          int indent, int level) -> std::forward_list<std::string> {
  BOOST_ASSERT(level >= 0);
  auto vertex_weights_str = dump_utils::dump_floating_point_array_as_braced_list( //
      base.coarsened.vertex_weights, indent, level + 1);
  return std::forward_list{
      fmt::format(".coarsened.adj_list = {}", //
                  dump_utils::dump_graph(base.coarsened.adj_list, indent, level + 1)),
      fmt::format(".coarsened.inv_adj_list = {}", //
                  dump_utils::dump_graph(base.coarsened.inv_adj_list, indent, level + 1)),
      fmt::format(".coarsened.vertex_weights = {}", std::move(vertex_weights_str)),
      fmt::format(".details = {}", base.details.dump(indent, level + 1)),
  };
}
} // namespace

auto WIMCoarsenedVertexBrief::dump(int indent, int level) const noexcept -> std::string {
  indent = std::max(indent, 0);
  level = std::max(level, 0);

  auto components = {
      fmt::format(".members = {}", members),
      fmt::format(".best_seed_index = {}", best_seed_index),
  };
  return dump_utils::merge_dumped_components(components, indent, level);
}

auto WBIMCoarsenedVertexBrief::dump(int indent, int level) const noexcept -> std::string {
  indent = std::max(indent, 0);
  level = std::max(level, 0);

  auto components = {
      fmt::format(".members = {}", members),
      fmt::format(".best_boosted_index = {}", best_boosted_index),
  };
  return dump_utils::merge_dumped_components(components, indent, level);
}

auto WIMCoarsenedVertexDetails::dump(int indent, int level) const noexcept -> std::string {
  indent = std::max(indent, 0);
  level = std::max(level, 0);

  auto components = coarsened_vertex_details_base_components(*this, indent, level);
  auto derived_components = std::forward_list{
      fmt::format(".heuristics_out_seed = {::.4f}", heuristics_out_seed),
      fmt::format(".p_seed_internal = {}", dump_p_internal(p_seed_internal, n_members(), indent, level + 1)),
      fmt::format(".best_seed_index = {}", best_seed_index),
  };
  components.splice_after(components.before_begin(), derived_components);
  components.sort(); // To put the members with the same prefix together
  return dump_utils::merge_dumped_components(components, indent, level);
}

auto WBIMCoarsenedVertexDetails::dump(int indent, int level) const noexcept -> std::string {
  indent = std::max(indent, 0);
  level = std::max(level, 0);

  auto components = coarsened_vertex_details_base_components(*this, indent, level);
  auto derived_components = std::forward_list{
      fmt::format(".heuristics_in_boost = {::.4f}", heuristics_in_boost),
      fmt::format(".p_boost_internal = {}", dump_p_internal(p_boost_internal, n_members(), indent, level + 1)),
      fmt::format(".best_boosted_index = {}", best_boosted_index),
      fmt::format(".is_seed = {}", is_seed),
  };
  components.splice_after(components.before_begin(), derived_components);
  components.sort(); // To put the members with the same prefix together
  return dump_utils::merge_dumped_components(components, indent, level);
}

auto WIMCoarsenedEdgeDetails::dump(int indent, int level) const noexcept -> std::string {
  indent = std::max(indent, 0);
  level = std::max(level, 0);

  auto p_seed_cross_str = dump_p_cross(p_seed_cross, n_members_left, n_members_right, indent, level + 1);
  auto components = coarsened_edge_details_base_components(*this, indent, level);
  auto derived_components = {
      fmt::format(".p_seed_cross = {}", std::move(p_seed_cross_str)),
      fmt::format(".merged = {}", merged),
  };
  components.splice_after(components.before_begin(), derived_components);
  components.sort(); // To put the members with the same prefix together
  return dump_utils::merge_dumped_components(components, indent, level);
}

auto WBIMCoarsenedEdgeDetails::dump(int indent, int level) const noexcept -> std::string {
  indent = std::max(indent, 0);
  level = std::max(level, 0);

  auto p_boost_cross_str = dump_p_cross(p_boost_cross, n_members_left, n_members_right, indent, level + 1);
  auto components = coarsened_edge_details_base_components(*this, indent, level);
  auto derived_components = std::forward_list{
      fmt::format(".p_boost_cross = {}", std::move(p_boost_cross_str)),
      fmt::format(".merged = {}", merged),
  };
  components.splice_after(components.before_begin(), derived_components);
  components.sort(); // To put the members with the same prefix together
  return dump_utils::merge_dumped_components(components, indent, level);
};

template <is_coarsened_vertex_brief VertexBriefType>
auto CoarseningBrief<VertexBriefType>::dump(int indent, int level) const noexcept -> std::string {
  indent = std::max(indent, 0);
  level = std::max(level, 0);

  auto components = coarsening_info_base_components(*this, indent, level);
  return dump_utils::merge_dumped_components(components, indent, level);
}

template <is_coarsened_vertex_details VertexDetailsType>
auto CoarseningDetails<VertexDetailsType>::dump(int indent, int level) const noexcept -> std::string {
  indent = std::max(indent, 0);
  level = std::max(level, 0);

  auto group_id_str = dump_utils::dump_integer_array_as_braced_list(group_id, indent, level + 1);
  auto index_in_group_str = dump_utils::dump_integer_array_as_braced_list(index_in_group, indent, level + 1);

  auto base_components = coarsening_info_base_components(*this, indent, level);
  auto components = std::forward_list{
      fmt::format(".group_id = {}", std::move(group_id_str)),
      fmt::format(".index_in_group = {}", std::move(index_in_group_str)),
  };
  components.splice_after(components.before_begin(), base_components);
  return dump_utils::merge_dumped_components(components, indent, level);
}

template <class BriefOrDetailsType>
auto WIMCoarsenGraphResult<BriefOrDetailsType>::dump(int indent, int level) const noexcept -> std::string {
  indent = std::max(indent, 0);
  level = std::max(level, 0);

  auto components = coarsen_graph_result_base_components(*this, indent, level);
  return dump_utils::merge_dumped_components(components, indent, level);
}

template <class BriefOrDetailsType>
auto WBIMCoarsenGraphResult<BriefOrDetailsType>::dump(int indent, int level) const noexcept -> std::string {
  indent = std::max(indent, 0);
  level = std::max(level, 0);

  auto base_components = coarsen_graph_result_base_components(*this, indent, level);
  auto components = std::forward_list{
      fmt::format(".coarsened_seeds = {}", //
                  dump_utils::dump_integer_array_as_braced_list(coarsened_seeds.vertex_list, indent, level + 1)),
  };
  components.splice_after(components.before_begin(), base_components);
  return dump_utils::merge_dumped_components(components, indent, level);
}

#define DEFINE_FREE_DUMP_FUNCTION_FOR_COARSENING_DETAILS(Type)                \
  auto dump(const Type& obj, int indent, int level) noexcept -> std::string { \
    return obj.dump(indent, level);                                           \
  }

COARSENING_DETAILS_TYPES(DEFINE_FREE_DUMP_FUNCTION_FOR_COARSENING_DETAILS)
#undef DECLARE_FREE_DUMP_FUNCTION_FOR_COARSENING_DETAILS
