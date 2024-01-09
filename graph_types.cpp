#include "graph_types.h"
#include <fmt/format.h>
#include <nameof.hpp>
#include <ylt/easylog.hpp>

namespace {
template <is_edge_property E>
constexpr auto serialization_header_chars() {
  constexpr auto type_name = NAMEOF_TYPE(E);
  // Aligns to multiple of 4: Given two positive integers x, y, ceil(x / y) = floor((x + y - 1) / y)
  // 2 : Two additional characters '[' and ']'. x = length + 2 above.
  // 3 : y = 4 above, then y - 1 = 3
  constexpr auto header_length = 4 * (type_name.length() + 2 + 3) / 4;
  auto res = std::array<char, header_length>{};
  res.front() = '[';
  auto [in, out] = ranges::copy(type_name, res.data() + 1);
  *out++ = ']';
  ranges::fill(out, res.end(), ' ');
  return res;
}

template <is_edge_property E>
auto read_directed_edge_list_generic(const std::string& input_file) noexcept -> rfl::Result<ReadGraphResult<E>> {
  auto fin = std::ifstream(input_file, std::ios::binary);
  if (!fin.is_open()) {
    return rfl::Error{fmt::format("Failed to open input file '{}'.", input_file)};
  }

  // Part 1: Type header
  constexpr auto expected_header = serialization_header_chars<E>();
  constexpr auto header_length = expected_header.size();
  ELOGFMT(DEBUG, "Expects header '{}' as type identifier of the graph to be deserialized.",
          std::string_view{expected_header.data(), header_length});

  auto header_buffer = std::array<char, header_length>{};
  fin.read(header_buffer.data(), header_length);
  if (header_buffer != expected_header) {
    constexpr auto msg_pattern = "Corrupted data: Header of edge type mismatch (expected: '{}').";
    return rfl::Error{fmt::format(msg_pattern, std::string_view{expected_header.data(), header_length})};
  }
  // Part 2: Size header: n = # of vertices
  vertex_id_t n;
  fin.read(reinterpret_cast<char*>(&n), sizeof(vertex_id_t));
  ELOGFMT(DEBUG, "In header of input file: n = {}", n);
  // Part 3: Vertex weights
  auto vertex_weights = std::vector<vertex_weight_t>(n);
  fin.read(reinterpret_cast<char*>(vertex_weights.data()), sizeof(vertex_weight_t) * n);
  // Part 4: Edges
  auto res = DirectedEdgeList<E>{};
  res.deserialize(fin);
  if (auto n_res = graph::num_vertices(res); n_res != n) {
    constexpr auto msg_pattern = "Corrupted data: # of vertices mismatch between the header (which is {}) "
                                 "and the graph body (which is {})";
    return rfl::Error{fmt::format(msg_pattern, n, n_res)};
  }

  return ReadGraphResult<E>{.edge_list = std::move(res), .vertex_weights = std::move(vertex_weights)};
}

template <is_edge_property E>
auto write_directed_edge_list_generic(const DirectedEdgeList<E>& graph, std::span<const vertex_weight_t> vertex_weights,
                                      const std::string& output_file) noexcept -> ResultVoid {
  auto fout = std::ofstream(output_file, std::ios::binary);
  if (!fout.is_open()) {
    return rfl::Error{fmt::format("Failed to open output file '{}'.", output_file)};
  }
  auto n = (vertex_id_t)graph::num_vertices(graph);
  if (vertex_weights.size() != n) {
    constexpr auto msg_pattern = "# of vertices mismatch between the graph (which is {}) "
                                 "and the vertex weight list (which is {})";
    return rfl::Error{fmt::format(msg_pattern, n, vertex_weights.size())};
  }

  constexpr auto header = serialization_header_chars<E>();
  ELOGFMT(DEBUG, "Using header '{}' as type identifier of serialized graph.",
          std::string_view{header.data(), header.size()});
  // Part 1: Type header
  fout.write(header.data(), header.size());
  // Part 2: Size header, n = # of vertices
  fout.write(reinterpret_cast<char*>(&n), sizeof(vertex_id_t));
  // Part 3: Vertex weights
  fout.write(reinterpret_cast<const char*>(vertex_weights.data()), sizeof(vertex_weight_t) * n);
  // Part 4: The graph edges
  graph.serialize(fout);
  return RESULT_VOID_SUCCESS;
}
} // namespace

auto read_directed_wim_edge_list(const std::string& input_file) noexcept -> rfl::Result<WIMReadGraphResult> {
  return read_directed_edge_list_generic<WIMEdge>(input_file);
}

auto read_directed_wbim_edge_list(const std::string& input_file) noexcept -> rfl::Result<WBIMReadGraphResult> {
  return read_directed_edge_list_generic<WBIMEdge>(input_file);
}

auto write_directed_wim_edge_list(const DirectedEdgeList<WIMEdge>& graph,
                                  std::span<const vertex_weight_t> vertex_weights,
                                  const std::string& output_file) noexcept -> ResultVoid {
  return write_directed_edge_list_generic(graph, vertex_weights, output_file);
}

auto write_directed_wbim_edge_list(const DirectedEdgeList<WBIMEdge>& graph,
                                   std::span<const vertex_weight_t> vertex_weights,
                                   const std::string& output_file) noexcept -> ResultVoid {
  return write_directed_edge_list_generic(graph, vertex_weights, output_file);
}
