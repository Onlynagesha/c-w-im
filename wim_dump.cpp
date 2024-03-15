#include "dump_utils.h"
#include "wim.h"
#include <fmt/ranges.h>
#include <magic_enum_format.hpp>

auto PRRSketch::dump(int indent, int level) const noexcept -> std::string {
  indent = std::max(indent, 0);
  level = std::max(level, 0);

  auto components = {
      fmt::format(".vertices = {}", //
                  dump_utils::dump_integer_array_as_braced_list(vertices, indent, level + 1)),
      fmt::format(".critical_vertices = {}", //
                  dump_utils::dump_integer_array_as_braced_list(critical_vertices, indent, level + 1)),
      fmt::format(".center = {}", center),
      fmt::format(".mapped_center = {}", mapped_center),
      fmt::format(".mapped_graph = {}", dump_utils::dump_graph(mapped_graph, indent, level + 1)),
      fmt::format(".inv_mapped_graph = {}", dump_utils::dump_graph(inv_mapped_graph, indent, level + 1)),
  };
  return dump_utils::merge_dumped_components(components, indent, level);
}

auto PRRSketchSet::dump(int indent, int level) const noexcept -> std::string {
  indent = std::max(indent, 0);
  level = std::max(level, 0);

  auto sketches_str = [&]() {
    auto to_str = TRANSFORM_VIEW(fmt::format("[{}] = {}", _1, sketches[_1].dump(indent, level + 2)));
    return dump_utils::merge_dumped_components( //
        range(sketches.size()) | to_str, indent, level + 1);
  }();
  auto inv_critical_str = [&]() {
    auto to_str = TRANSFORM_VIEW(fmt::format("{}", _1));
    return dump_utils::merge_dumped_components(inv_critical | to_str, indent, level + 1);
  }();

  auto components = {
      fmt::format(".graph = {:p}, |V| = {}, |E| = {}", //
                  static_cast<const void*>(graph), graph::num_vertices(*graph), graph->num_edges()),
      fmt::format(".inv_graph = {:p}, |V| = {}, |E| = {}", //
                  static_cast<const void*>(inv_graph), graph::num_vertices(*inv_graph), graph->num_edges()),
      fmt::format(".seeds = {}", //
                  dump_utils::dump_integer_array_as_braced_list(seeds->vertex_list, indent, level + 1)),
      fmt::format(".sketches = {}", sketches_str),
      fmt::format(".inv_critical = {}", inv_critical_str),
  };
  return dump_utils::merge_dumped_components(components, indent, level);
}

#define DEFINE_FREE_DUMP_FUNCTION_FOR_WIM_TYPES(Type)                         \
  auto dump(const Type& obj, int indent, int level) noexcept -> std::string { \
    return obj.dump(indent, level);                                           \
  }

WIM_DUMP_TYPES(DEFINE_FREE_DUMP_FUNCTION_FOR_WIM_TYPES)
#undef DEFINE_FREE_DUMP_FUNCTION_FOR_WIM_TYPES
