#pragma once

#include "experiments.h"
#include "graph_connectivity.h"
#include "utils/easylog.h"
#include <nlohmann/json.hpp>

namespace exp_io {
inline auto init_easylog(const CommonExperimentParams& params) {
  easylog::set_min_severity(params.log_severity);
  // Note: async=true may trigger a bug that the program fails to terminate after everything is finished.
  if (!params.log_output_file.empty()) {
    easylog::init_log(params.log_severity, params.log_output_file, false, params.log_console);
  }
}

// Returns nullptr if out_path is empty
inline auto create_fout_ptr(const std::string& out_path) -> std::unique_ptr<std::ofstream> {
  auto fout_ptr = std::unique_ptr<std::ofstream>{};
  if (!out_path.empty()) {
    auto fout = std::ofstream{out_path};
    if (!fout.is_open()) {
      ELOGFMT(ERROR, "Failed to open output file '{}'. Uses log output as fallback.", out_path);
    } else {
      ELOGFMT(INFO, "Successfully opens output file '{}'.", out_path);
      fout_ptr = std::make_unique<std::ofstream>(std::move(fout));
    }
  } else {
    ELOG_INFO << "No output file specified. Uses log output as fallback.";
  }
  return fout_ptr;
}

inline auto dump_to_fout_ptr(std::ofstream* fout, std::string_view contents) -> void {
  if (fout != nullptr) {
    MYLOG_FMT_DEBUG("JSON output: {}", contents);
    (*fout) << contents;
  } else {
    ELOGFMT(INFO, "JSON output: {}", contents);
  }
}

template <is_edge_property E>
inline auto dump_vertices_selected(const AdjacencyListPair<E>& graph, std::span<const vertex_id_t> vertices_selected)
    -> std::string {
  auto [n, m] = graph.graph_n_m();
  auto res = fmt::format("Vertices selected, with |E| / |V| = {:.3f}:", 1.0 * m / n);
  for (auto v : vertices_selected) {
    BOOST_ASSERT_MSG(v < n, "Vertex index out of range [0, n).");
    res += fmt::format("\n\tvertex-id = {}, in-degree = {}, out-degree = {}", //
                       v, graph.in_degree(v), graph.out_degree(v));
  }
  return res;
}

template <is_edge_property E>
inline auto write_graph_basic_information(json& json_root, const AdjacencyList<E>& graph) {
  auto n = graph.num_vertices()[0];
  auto m = graph.num_edges();
  json_root["n"] = n;
  json_root["m"] = m;
  return std::tuple{n, m};
}

struct GraphConnectivityInfo {
  vertex_id_t n_wcc;
  vertex_id_t n_scc;
};

template <is_edge_property E>
inline auto write_graph_connectivity_information(json& json_root, const AdjacencyListPair<E>& graph)
    -> GraphConnectivityInfo {
  auto res = GraphConnectivityInfo{
      .n_wcc = n_weakly_connected_components(graph.adj_list, graph.inv_adj_list),
      .n_scc = n_strongly_connected_components(graph.adj_list),
  };
  json_root["n_wcc"] = res.n_wcc;
  json_root["n_scc"] = res.n_scc;
  return res;
}
} // namespace exp_io
