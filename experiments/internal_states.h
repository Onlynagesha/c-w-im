#pragma once

#include "coarsening_details_types.h"
#include "graph_types.h"
#include <nlohmann/json.hpp>

#define EXP_STATES_TO_JSON_FUNCTION_AUTO_GENERATED(F) \
  F(WIMSketchingInfo)                                 \
  F(WIMSketchingGetSeedsResult)                       \
  F(WBIMSketchingInfo)                                \
  F(WBIMSketchingGetBoostedResult)                    \
  F(WIMSimulationInfo)                                \
  F(WBIMSimulationInfo)                               \
  F(WIMSketchingSimulationResult)                     \
  F(WBIMSketchingSimulationResult)                    \
  F(CoarseningInfo)                                   \
  F(WIMExpansionInfo)                                 \
  F(WBIMExpansionInfo)                                \
  F(WIMExpansionResult)                               \
  F(WBIMExpansionResult)

namespace exp_states {
using VertexList = std::vector<vertex_id_t>;
using VertexListList = std::vector<VertexList>;

template <class T>
  requires(requires(const T& obj) {
    { to_json(obj) } -> std::same_as<json>;
  })
inline auto to_json(const std::vector<T>& vector) -> json {
  auto res = json{};
  for (const auto& obj : vector) {
    res.push_back(to_json(obj));
  }
  return res;
}

struct WIMSketchingInfo {
  size_t n_sketches;
  VertexList selected_seeds;

  double average_sketch_size;
  double sketching_total_time_usage;
  double seed_selecting_time_usage;
};

struct WIMSketchingGetSeedsResult {
  VertexListList selected_seeds;
  std::vector<WIMSketchingInfo> sketching_info;
};

struct WBIMSketchingInfo {
  size_t n_sketches;
  VertexList selected_boosted;
  VertexList selected_boosted_by_critical;

  double average_sketch_n_vertices;
  double average_sketch_n_edges;
  size_t total_sketch_size_bytes;
  double sketching_success_rate;

  double sketching_total_time_usage;
  double boosted_selecting_time_usage;
  double boosted_selecting_time_usage_by_critical;
};

struct WBIMSketchingGetBoostedResult {
  VertexListList selected_boosted;
  VertexListList selected_boosted_by_critical;
  std::vector<WBIMSketchingInfo> sketching_info;
};

struct WIMSimulationInfo {
  vertex_id_t n_seeds;
  double objective_function;
  double time_usage;
};

struct WBIMSimulationInfo {
  vertex_id_t n_boosted;
  double objective_function; // Takes the difference: F(B;S) - F(EmptySet;S)
  double time_usage;
};

struct WIMSketchingSimulationResult {
  uint64_t n_sketches;
  std::vector<WIMSimulationInfo> simulation_results;
};

struct WBIMSketchingSimulationResult {
  uint64_t n_sketches;
  std::vector<WBIMSimulationInfo> simulation_results;
};

struct CoarseningInfo {
  vertex_id_t n_coarsened;
  vertex_id_t m_coarsened;
  double time_usage;
};

template <is_edge_property E>
struct CoarseningResult {
  using BriefResultType = typename CoarseningDataTraits<E>::BriefResult;

  std::vector<BriefResultType> coarsen_results;
  CoarseningInfo info;
};

template <is_edge_property E>
inline auto exp_state_to_json(const CoarseningResult<E>& obj) -> json {
  return exp_state_to_json(obj.info);
}

struct WIMExpansionInfo {
  VertexList expanded_seeds;
  double time_usage;
};

struct WBIMExpansionInfo {
  VertexList expanded_boosted;
  double time_usage;
};

struct WIMExpansionResult {
  VertexListList expanded_seeds;
  std::vector<double> time_usage;

  auto info_item(size_t index) const -> WIMExpansionInfo {
    BOOST_ASSERT(index < expanded_seeds.size() && index < time_usage.size());
    return WIMExpansionInfo{
        .expanded_seeds = expanded_seeds[index],
        .time_usage = time_usage[index],
    };
  }
};

struct WBIMExpansionResult {
  VertexListList expanded_boosted;
  std::vector<double> time_usage;

  auto info_item(size_t index) const -> WBIMExpansionInfo {
    BOOST_ASSERT(index < expanded_boosted.size() && index < time_usage.size());
    return WBIMExpansionInfo{
        .expanded_boosted = expanded_boosted[index],
        .time_usage = time_usage[index],
    };
  }
};

template <is_edge_property E>
struct ExpansionTraits {
  static_assert(rfl::always_false_v<E>, "Invalid edge type.");
};

template <>
struct ExpansionTraits<WIMEdge> {
  using InfoType = WIMExpansionInfo;
  using ResultType = WIMExpansionResult;
};

template <>
struct ExpansionTraits<WBIMEdge> {
  using InfoType = WBIMExpansionInfo;
  using ResultType = WBIMExpansionResult;
};

#define DECLARE_EXP_STATE_TO_JSON_FUNCTION(Type) auto to_json(const Type& obj) noexcept -> json;
// Functions declares in namespace exp_states
EXP_STATES_TO_JSON_FUNCTION_AUTO_GENERATED(DECLARE_EXP_STATE_TO_JSON_FUNCTION)
} // namespace exp_states

#undef DECLARE_EXP_STATE_TO_JSON_FUNCTION
