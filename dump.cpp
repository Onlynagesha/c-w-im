#include "dump.h"
#include "create_dataset.h"
#include "graph_types.h"
#include "utils/demangle.h"
#include <fmt/ranges.h>
#include <nlohmann/json.hpp>
#include <nwgraph/adaptors/edge_range.hpp>
#include <nwgraph/util/demangle.hpp>
#include <rfl.hpp>
#include <rfl/json.hpp>

namespace {
template <class T>
inline auto dump_generic(const T& value, int indent = 0) noexcept -> std::string try {
  auto unindented = rfl::json::write(value);
  if (indent <= 0) {
    return unindented;
  }
  auto json_obj = json::parse(unindented);
  return json_obj.dump(indent);
} catch (std::exception& e) {
  constexpr auto msg_pattern = "<DUMP ERROR: [{}] {}>";
  return fmt::format(msg_pattern, demangle_type_name(e), e.what());
} catch (...) {
  return "<DUMP ERROR: Unknown>";
}

template <class T>
inline auto dump_as_json_generic(const T& value) noexcept -> json {
  auto text = rfl::json::write(value);
  return json::parse(text);
}

auto merge_dumped_components(std::span<const std::string> components, int outer_indent, int inner_indent)
    -> std::string {
  if (inner_indent <= 0) {
    return fmt::format("{{{}}}", components);
  }
  auto res = "{"s;
  for (const auto& c : components) {
    res += '\n' + std::string(outer_indent + inner_indent, ' ');
    res += c;
  }
  res += '\n' + std::string(outer_indent, ' ') + '}';
  return res;
}

auto dump_impl(const CoarsenedVertexDetails& details, int outer_indent, int inner_indent) noexcept -> std::string {
  outer_indent = std::max(outer_indent, 0);
  inner_indent = std::max(inner_indent, 0);
  using PInternalContainer = CoarsenedVertexDetails::PInternalContainer;
  auto dump_p_internal = [&](std::string_view name, const PInternalContainer& p_internal) {
    auto m = details.members.size();
    auto rows = range(m) | views::transform(LAMBDA_1(fmt::format("{::.4f}", p_internal[_1] | views::take(m))));
    if (inner_indent <= 0) {
      return fmt::format(".{} = {}", name, rows);
    }
    auto res = fmt::format(".{} = [", name);
    for (auto i : range(m)) {
      res += '\n' + std::string(outer_indent + inner_indent * 2, ' ');
      res += rows[i];
    }
    res += '\n' + std::string(outer_indent + inner_indent, ' ') + ']';
    return res;
  };
  auto components = {fmt::format(".members = {}", details.members),
                     fmt::format(".vertex_weights = {::.4f}", details.vertex_weights),
                     fmt::format(".heuristics_in = {::.4f}", details.heuristics_in),
                     fmt::format(".heuristics_out = {::.4f}", details.heuristics_out),
                     fmt::format(".heuristics_out_seed = {::.4f}", details.heuristics_out_seed),
                     dump_p_internal("p_internal", details.p_internal),
                     dump_p_internal("p_seed_internal", details.p_seed_internal)};
  return merge_dumped_components(components, outer_indent, inner_indent);
}

auto dump_impl(const CoarsenedEdgeDetails& details, int outer_indent, int inner_indent) noexcept -> std::string {
  outer_indent = std::max(outer_indent, 0);
  inner_indent = std::max(inner_indent, 0);
  using PCrossContainer = CoarsenedEdgeDetails::PCrossContainer;
  auto dump_p_cross = [&](std::string_view name, const PCrossContainer& p_cross) {
    auto n_rows = details.n_members_left;
    auto n_cols = details.n_members_right;
    auto rows = range(n_rows) | views::transform(LAMBDA_1(fmt::format("{::.4f}", p_cross[_1] | views::take(n_cols))));
    if (inner_indent <= 0) {
      return fmt::format(".{} = {}", name, rows);
    }
    auto res = fmt::format(".{} = [", name);
    for (auto i : range(n_rows)) {
      res += '\n' + std::string(outer_indent + inner_indent * 2, ' ');
      res += rows[i];
    }
    res += '\n' + std::string(outer_indent + inner_indent, ' ') + ']';
    return res;
  };
  auto components = {dump_p_cross("p_cross", details.p_cross), dump_p_cross("p_seed_cross", details.p_seed_cross),
                     fmt::format(".merged = {}", details.merged)};
  return merge_dumped_components(components, outer_indent, inner_indent);
}

auto dump_impl(CoarsenGraphResult& result, int indent) noexcept -> std::string {
  auto res = "coarsened_graph = ["s;
  auto first = true;
  for (auto [u, v, w] : graph::make_edge_range<0>(result.coarsened_graph)) {
    if (!first) {
      res += ',';
    } else {
      first = false;
    }
    if (indent > 0) {
      res += '\n' + std::string(indent, ' ');
    }
    res += fmt::format("{} -> {}: {{.p = {:.4f}, .p_seed = {:.4f}}}", u, v, w.p, w.p_seed);
  }
  res += fmt::format("\n],\ncoarsened_vertex_weights = {::.4f}", result.coarsened_vertex_weights);
  res += "\ndetails:\n" + dump(result.details, indent);
  return res;
}
} // namespace

#define IMPLEMENT_DUMP_FUNCTIONS(Type)                               \
  auto dump(const Type& value, int indent) noexcept -> std::string { \
    return ::dump_generic(value, indent);                            \
  }                                                                  \
  auto dump_as_json(const Type& value) noexcept -> json {            \
    return ::dump_as_json_generic(value);                            \
  }

DUMP_REGISTERED_TYPES_AUTO_GENERATED(IMPLEMENT_DUMP_FUNCTIONS)
#undef IMPLEMENT_DUMP_FUNCTIONS

// ---- Manual implementation ----

auto dump(const CoarsenedVertexDetails& details, int indent) noexcept -> std::string {
  return dump_impl(details, indent, indent);
}

auto dump(const CoarsenedEdgeDetails& details, int indent) noexcept -> std::string {
  return dump_impl(details, indent, indent);
}

auto dump(const CoarseningDetails& details, int indent) noexcept -> std::string {
  auto res = fmt::format("n = {}\nn_groups = {}", details.n, details.n_groups);
  res += '\n' + dump_array("group_id", details.group_id);
  res += '\n' + dump_array("index_in_group", details.index_in_group);
  // groups
  res += "\ngroups = [";
  for (auto [i, g] : details.groups | views::enumerate) {
    if (i != 0) {
      res += ',';
    }
    res += fmt::format("\n{}Group #{}: {}", std::string(indent, ' '), i, dump_impl(g, indent, indent));
  }
  res += "\n],";
  // edges
  res += "\nedges = [";
  bool first = true;
  // for (const auto& [p, e] : details.edges) {
  //   if (!first) {
  //     res += ',';
  //   } else {
  //     first = false;
  //   }
  //   res += fmt::format("\n{}Edge {}: {}", std::string(indent, ' '), p, dump_impl(e, indent, indent));
  // }
  res += "\n]";
  return res;
}

auto dump(const CoarsenGraphResult& result, int indent) noexcept -> std::string {
  // Workaround: const reference triggers a bug in nwgraph. The graph is not changed actually.
  return dump_impl(const_cast<CoarsenGraphResult&>(result), indent);
}
