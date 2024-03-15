#pragma once

#include "graph_types.h"
#include "utils/result.h"
#include "utils/static_vector.h"
#include "wim.h"
#include <array>
#include <rfl/Rename.hpp>
#include <rfl/Validator.hpp>
#include <rfl/comparisons.hpp>

#define COARSENING_DETAILS_TYPES(F) \
  F(WIMCoarsenedVertexBrief)        \
  F(WIMCoarsenedVertexDetails)      \
  F(WIMCoarsenedEdgeDetails)        \
  F(WIMCoarseningBrief)             \
  F(WIMCoarseningDetails)           \
  F(WIMCoarsenGraphBriefResult)     \
  F(WIMCoarsenGraphDetailedResult)

#define DECLARE_FREE_DUMP_FUNCTION_FOR_COARSENING_DETAILS(Type) \
  auto dump(const Type& obj, int indent = 0, int level = 0) noexcept -> std::string;

enum class NeighborMatchRule { HEM_P_MAX, HEM_P_PRODUCT, LEM_P_MAX, LEM_P_PRODUCT };

enum class EdgeWeightRule { SEPARATE_SIMPLE, MERGED_SIMPLE, SEPARATE_PRECISE, MERGED_PRECISE };

enum class EdgeSeedWeightRule { AVERAGE, MAX, BEST_SEED_INDEX };

enum class InOutHeuristicRule { UNIT, COUNT, P, W, SEED };

enum class VertexWeightRule { AVERAGE, AVERAGE_BY_PATHS, SUM };

enum class SeedExpandingRule { LOCAL, SIMULATIVE, ITERATIVE };

struct CoarseningParams {
  NeighborMatchRule neighbor_match_rule = NeighborMatchRule::HEM_P_MAX;
  EdgeWeightRule edge_weight_rule = EdgeWeightRule::SEPARATE_PRECISE;
  EdgeSeedWeightRule edge_seed_weight_rule = EdgeSeedWeightRule::BEST_SEED_INDEX;
  InOutHeuristicRule in_out_heuristic_rule = InOutHeuristicRule::P;
  VertexWeightRule vertex_weight_rule = VertexWeightRule::AVERAGE_BY_PATHS;
  vertex_id_t max_distance_from_seed = 6;
};

struct ExpandingParams {
  SeedExpandingRule seed_expanding_rule;
  // Used for MERGED only
  using NIterations = rfl::Rename<              //
      "expanding_n_iterations",                 //
      rfl::Validator<uint64_t, rfl::Minimum<1>> //
      >;
  NIterations n_iterations = 10;

  using SimulationTryCount = rfl::Rename<       //
      "expanding_simulation_try_count",         //
      rfl::Validator<uint64_t, rfl::Minimum<1>> //
      >;
  SimulationTryCount simulation_try_count = 1'000;
};

// ---- Step 1: Converts to undirected edge ----

auto merge_wim_edge_to_undirected(const AdjacencyList<WIMEdge>& graph, const CoarseningParams& params) noexcept
    -> AdjacencyList<edge_probability_t>;

auto merge_wbim_edge_to_undirected(const AdjacencyList<WBIMEdge>& graph, const CoarseningParams& params) noexcept
    -> AdjacencyList<edge_probability_t>;

template <is_edge_property E>
inline auto merge_edge_to_undirected(const AdjacencyList<E>& graph, const CoarseningParams& params) noexcept {
  if constexpr (std::is_same_v<E, WIMEdge>) {
    return merge_wim_edge_to_undirected(graph, params);
  } else if constexpr (std::is_same_v<E, WBIMEdge>) {
    return merge_wbim_edge_to_undirected(graph, params);
  } else {
    static_assert(rfl::always_false_v<E>, "Invalid edge type.");
  }
}

// ---- Step 2: Performs Mongoose matching algorithm ----

struct MongooseMatchResult {
  vertex_id_t n_groups;
  std::vector<vertex_id_t> group_id;
};

auto mongoose_match(const AdjacencyList<edge_probability_t>& bidir_graph, const CoarseningParams& params) noexcept
    -> MongooseMatchResult;

auto mongoose_match_s(const AdjacencyList<edge_probability_t>& bidir_graph, std::span<const vertex_id_t> seeds,
                      const CoarseningParams& params) noexcept -> MongooseMatchResult;

// ---- Step 3: Performs coarsening by the groups obtained previously ----

struct WIMCoarsenedVertexBrief {
  // Mongoose algorithm ensures that the size of each group is at most 3 (match + 1 adopted)
  static constexpr auto MAX_N_MEMBERS = size_t{3};
  using MemberContainer = StaticVector<vertex_id_t, MAX_N_MEMBERS>;

  // List of size M, # of members in the groups, i.e. the vertices before coarsening into current group
  MemberContainer members;
  // Index in range [0, M), the best expanded seed by local information.
  size_t best_seed_index;

  auto dump(int indent = 0, int level = 0) const noexcept -> std::string;
};

struct WIMCoarsenedVertexDetails {
  static constexpr auto MAX_N_MEMBERS = WIMCoarsenedVertexBrief::MAX_N_MEMBERS;

  using MemberContainer = WIMCoarsenedVertexBrief::MemberContainer;
  using VertexWeightsContainer = StaticVector<vertex_weight_t, MAX_N_MEMBERS>;
  using HeuristicsContainer = StaticVector<edge_probability_t, MAX_N_MEMBERS>;
  using PInternalContainer = Array2D<edge_probability_t, MAX_N_MEMBERS, MAX_N_MEMBERS>;

  // List of size M, # of members in the groups, i.e. the vertices before coarsening into current group
  MemberContainer members;
  // List of size M. vertex_weights[i] = Vertex weight of members[i]
  VertexWeightsContainer vertex_weights;
  // List of size M. heuristics_in[i] = In-heuristic value of members[i]
  // Note: NOT NORMALIZED; MAY BE ALL-ZERO
  HeuristicsContainer heuristics_in;
  // List of size M. heuristics_out[i] = Out-heuristic value of members[i]
  // Note: NOT NORMALIZED; MAY BE ALL-ZERO
  HeuristicsContainer heuristics_out;
  // List of size M. heuristics_out[i] = Out-heuristic value of members[i] as seed
  // Note: NOT NORMALIZED; MAY BE ALL-ZERO
  HeuristicsContainer heuristics_out_seed;
  // Matrix of shape (M, M). p_internal[i][j] = p of the directed edge members[i] -> members[j]
  PInternalContainer p_internal;
  // Matrix of shape (M, M). p_seed_internal[i][j] = p_seed of the directed edge members[i] -> members[j]
  PInternalContainer p_seed_internal;
  // Index in range [0, M), the best expanded seed by local information.
  size_t best_seed_index;

  auto dump(int indent = 0, int level = 0) const noexcept -> std::string;

  auto n_members() const -> size_t {
    return members.size();
  }

  auto member_indices() const -> ranges::iota_view<size_t, size_t> {
    return range(n_members());
  }

  template <size_t N>
  auto get_members() const -> std::array<vertex_id_t, N> {
    BOOST_ASSERT_MSG(members.size() == N, "# of members mismatch.");
    auto res = std::array<vertex_id_t, N>{}; // Supports structured binding
    ranges::copy(members, res.begin());
    return res;
  }
};

/*
Details of the coarsened edge Gu -> Gv,
where the coarsened vertex Gu -> {u1, u2, ...}, Gv -> {v1, v2, ...} after expansion.
*/
struct WIMCoarsenedEdgeDetails {
  static constexpr auto MAX_N_MEMBERS = WIMCoarsenedVertexDetails::MAX_N_MEMBERS;
  using PCrossContainer = Array2D<edge_probability_t, MAX_N_MEMBERS, MAX_N_MEMBERS>;

  // Size of {u1, u2, ...}. Demoted as Mu
  size_t n_members_left;
  // Size of {v1, v2, ...}. Denoted as Mv
  size_t n_members_right;
  // p_cross[i][j] = p(ui, vj)
  PCrossContainer p_cross;
  // p_seed_cross[i][j] = p_seed(ui, vj)
  PCrossContainer p_seed_cross;
  // Properties of the merged edge
  WIMEdge merged;

  auto dump(int indent = 0, int level = 0) const noexcept -> std::string;
};

struct WIMCoarseningBrief {
  // N, # of vertices before coarsening
  vertex_id_t n;
  // Nc, # of vertices in the coarsened graph.
  vertex_id_t n_coarsened;
  // List of size Nc. groups[g] represents the group {v1, v2, ...} whose group index is g,
  // i.e. all the vertices v1, v2, ... that are coarsened to g.
  std::vector<WIMCoarsenedVertexBrief> groups;

  auto dump(int indent = 0, int level = 0) const noexcept -> std::string;
};

struct WIMCoarseningDetails {
  // N, # of vertices before coarsening
  vertex_id_t n;
  // Nc, # of vertices in the coarsened graph.
  vertex_id_t n_coarsened;
  // List of size N.
  std::vector<vertex_id_t> group_id;
  // List of size N.
  // e.g. Let {v1, v2, v3} be a group, with group_id[v1] = group_id[v2] = group_id[v3] = g, v1 < v2 < v3,
  // Then index_in_group[v1] = 0, index_in_group[v2] = 1, index_in_group[v3] = 2
  std::vector<vertex_id_t> index_in_group;
  // List of size Nc. groups[g] represents the group {v1, v2, v3} whose group index is g.
  std::vector<WIMCoarsenedVertexDetails> groups;

  auto dump(int indent = 0, int level = 0) const noexcept -> std::string;

private:
  template <class GetFn, std::convertible_to<vertex_id_t>... Args>
  auto to_get_fn_result_generic(GetFn&& get_fn, Args... vertex_indices) const {
    if constexpr (sizeof...(Args) == 1) {
      return get_fn(vertex_indices...);
    } else {
      return std::tuple{get_fn(vertex_indices)...};
    }
  }

public:
  template <std::convertible_to<vertex_id_t>... Args>
  auto to_group_id(Args... vertex_indices) const {
    BOOST_ASSERT_MSG(((vertex_indices < n) && ...), "Vertex indices are out of range [0, n).");
    auto get_fn = LAMBDA_1(group_id[_1]);
    return to_get_fn_result_generic(get_fn, vertex_indices...);
  }

  template <std::convertible_to<vertex_id_t>... Args>
  auto to_group_ptr(Args... vertex_indices) const {
    BOOST_ASSERT_MSG(((vertex_indices < n) && ...), "Vertex indices are out of range [0, n).");
    auto get_fn = LAMBDA_1(groups.data() + group_id[_1]);
    return to_get_fn_result_generic(get_fn, vertex_indices...);
  }

  template <std::convertible_to<vertex_id_t>... Args>
  auto to_index_in_group(Args... vertex_indices) const {
    BOOST_ASSERT_MSG(((vertex_indices < n) && ...), "Vertex indices are out of range [0, n).");
    auto get_fn = LAMBDA_1(index_in_group[_1]);
    return to_get_fn_result_generic(get_fn, vertex_indices...);
  }

  template <std::convertible_to<vertex_id_t>... Args>
  auto to_group_size(Args... vertex_indices) const {
    BOOST_ASSERT_MSG(((vertex_indices < n) && ...), "Vertex indices are out of range [0, n).");
    auto get_fn = LAMBDA_1(groups[group_id[_1]].n_members());
    return to_get_fn_result_generic(get_fn, vertex_indices...);
  }

  template <std::convertible_to<vertex_id_t>... Args>
  auto group_id_to_size(Args... group_indices) const {
    BOOST_ASSERT_MSG(((group_indices < n_coarsened) && ...), "Group indices are out of range [0, nc).");
    auto get_fn = LAMBDA_1(groups[_1].n_members());
    return to_get_fn_result_generic(get_fn, group_indices...);
  }
};

template <same_as_either<WIMCoarseningBrief, WIMCoarseningDetails> DetailsType>
struct WIMCoarsenGraphResult {
  WIMAdjacencyListPair coarsened;
  DetailsType details;

  auto dump(int indent = 0, int level = 0) const noexcept -> std::string;
};

using WIMCoarsenGraphBriefResult = WIMCoarsenGraphResult<WIMCoarseningBrief>;
using WIMCoarsenGraphDetailedResult = WIMCoarsenGraphResult<WIMCoarseningDetails>;

auto coarsen_wim_graph_by_match(const AdjacencyList<WIMEdge>& graph, const InvAdjacencyList<WIMEdge>& inv_graph,
                                std::span<const vertex_weight_t> vertex_weights, vertex_id_t n_groups,
                                std::span<const vertex_id_t> group_id, const CoarseningParams& params) noexcept
    -> rfl::Result<WIMCoarsenGraphBriefResult>;

auto coarsen_wim_graph_by_match_d(const AdjacencyList<WIMEdge>& graph, const InvAdjacencyList<WIMEdge>& inv_graph,
                                  std::span<const vertex_weight_t> vertex_weights, vertex_id_t n_groups,
                                  std::span<const vertex_id_t> group_id, const CoarseningParams& params) noexcept
    -> rfl::Result<WIMCoarsenGraphDetailedResult>;

// ---- Combination of Step 1 - 3 ----

auto coarsen_wim_graph(const AdjacencyList<WIMEdge>& graph, const InvAdjacencyList<WIMEdge>& inv_graph,
                       std::span<const vertex_weight_t> vertex_weights, const CoarseningParams& params) noexcept
    -> rfl::Result<WIMCoarsenGraphBriefResult>;

inline auto coarsen_wim_graph_p(const WIMAdjacencyListPair& graph, const CoarseningParams& params) noexcept
    -> rfl::Result<WIMCoarsenGraphBriefResult> {
  return coarsen_wim_graph(graph.adj_list, graph.inv_adj_list, graph.vertex_weights, params);
}

auto coarsen_wim_graph_d(const AdjacencyList<WIMEdge>& graph, const InvAdjacencyList<WIMEdge>& inv_graph,
                         std::span<const vertex_weight_t> vertex_weights, const CoarseningParams& params) noexcept
    -> rfl::Result<WIMCoarsenGraphDetailedResult>;

inline auto coarsen_wim_graph_d_p(const WIMAdjacencyListPair& graph, const CoarseningParams& params) noexcept
    -> rfl::Result<WIMCoarsenGraphDetailedResult> {
  return coarsen_wim_graph_d(graph.adj_list, graph.inv_adj_list, graph.vertex_weights, params);
}

// ---- Step 4: Expanding seeds ----

struct SelectBestSeedResult {
  size_t index_in_group;
};

auto select_best_seed_in_group(const WIMCoarsenedVertexDetails& v, InOutHeuristicRule rule) noexcept
    -> SelectBestSeedResult;

struct ExpandSeedResult {
  std::vector<vertex_id_t> expanded_seeds;
};

auto expand_wim_seed_vertices(const AdjacencyList<WIMEdge>& graph, std::span<const vertex_weight_t> vertex_weights,
                              const WIMCoarsenGraphBriefResult& coarsening_result,
                              std::span<const vertex_id_t> coarsened_seeds, const ExpandingParams& params) noexcept
    -> rfl::Result<ExpandSeedResult>;

auto expand_wim_seed_vertices_d(const AdjacencyList<WIMEdge>& graph, std::span<const vertex_weight_t> vertex_weights,
                                const WIMCoarsenGraphDetailedResult& coarsening_result,
                                std::span<const vertex_id_t> coarsened_seeds, const ExpandingParams& params) noexcept
    -> rfl::Result<ExpandSeedResult>;

inline auto further_expand_wim_seed_vertices(const WIMCoarsenGraphBriefResult& last_result,
                                             const WIMCoarsenGraphBriefResult& cur_result,
                                             std::span<const vertex_id_t> coarsened_seeds,
                                             const ExpandingParams& params) -> rfl::Result<ExpandSeedResult> {
  return expand_wim_seed_vertices( //
      last_result.coarsened.adj_list, last_result.coarsened.vertex_weights, cur_result, coarsened_seeds, params);
}

inline auto further_expand_wim_seed_vertices_d(const WIMCoarsenGraphDetailedResult& last_result,
                                               const WIMCoarsenGraphDetailedResult& cur_result,
                                               std::span<const vertex_id_t> coarsened_seeds,
                                               const ExpandingParams& params) -> rfl::Result<ExpandSeedResult> {
  return expand_wim_seed_vertices_d( //
      last_result.coarsened.adj_list, last_result.coarsened.vertex_weights, cur_result, coarsened_seeds, params);
}

// ---- Free functions for dumping ----

COARSENING_DETAILS_TYPES(DECLARE_FREE_DUMP_FUNCTION_FOR_COARSENING_DETAILS)

#undef DECLARE_FREE_DUMP_FUNCTION_FOR_COARSENING_DETAILS
