#pragma once

#include "graph_types.h"
#include "utils/result.h"
#include "utils/static_vector.h"
#include <array>
#include <rfl/Rename.hpp>
#include <rfl/Validator.hpp>
#include <rfl/comparisons.hpp>

enum class NeighborMatchRule { HEM_P_MAX, HEM_P_PRODUCT, LEM_P_MAX, LEM_P_PRODUCT };

enum class GroupPathProbabilityRule { P_SEPARATE, P_MERGED };

enum class GroupInOutRule { NONE, P, W };

enum class VertexWeightRule { AVERAGE, AVERAGE_BY_PATHS, SUM };
// C-w-BIM only
enum class SeedMergeRule { UNUSED, S_SINGLE, S_MERGED };

enum class SeedExpansionRule { S_LOCAL, S_SIMULATIVE, S_ITERATIVE };

struct CoarseningParams {
  NeighborMatchRule neighbor_match_rule;
  GroupPathProbabilityRule group_path_probability_rule;
  GroupInOutRule group_in_out_rule;
  VertexWeightRule vertex_weight_rule;
  SeedMergeRule seed_merge_rule;
};

struct ExpandingParams {
  SeedExpansionRule seed_expansion_rule;
  // Unused for S_SEPARATE only
  using NIterations = rfl::Rename<"expanding_n_iterations", rfl::Validator<uint64_t, rfl::Minimum<1>>>;
  NIterations n_iterations = 10;

  using SimulationTryCount = rfl::Rename<"expanding_simulation_try_count", rfl::Validator<uint64_t, rfl::Minimum<1>>>;
  SimulationTryCount simulation_try_count = 10'000;
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

auto mongoose_match(const AdjacencyList<edge_probability_t>& graph, const CoarseningParams& params) noexcept
    -> MongooseMatchResult;

auto mongoose_match_with_seeds(const AdjacencyList<edge_probability_t>& graph, std::span<const vertex_id_t> seeds,
                               const CoarseningParams& params) noexcept -> MongooseMatchResult;

// ---- Step 3: Performs coarsening by the groups obtained previously ----

struct CoarsenedVertexDetails {
  // Mongoose algorithm ensures that the size of each group is at most 3 (match + 1 adopted)
  static constexpr auto MAX_NUM_MEMBERS = size_t{3};
  using MemberContainer = StaticVector<vertex_id_t, MAX_NUM_MEMBERS>;
  using VertexWeightsContainer = StaticVector<vertex_weight_t, MAX_NUM_MEMBERS>;
  using HeuristicsContainer = StaticVector<edge_probability_t, MAX_NUM_MEMBERS>;
  using PInternalContainer = Array2D<edge_probability_t, MAX_NUM_MEMBERS, MAX_NUM_MEMBERS>;

  // List of size m, # of members in the groups, i.e. the vertices before coarsening into current group
  MemberContainer members;
  // List of size m. vertex_weights[i] = Vertex weight of members[i]
  VertexWeightsContainer vertex_weights;
  // List of size m. heuristics_in[i] = In-heuristic value of members[i]
  // Note: NOT NORMALIZED; MAY BE ALL-ZERO
  HeuristicsContainer heuristics_in;
  // List of size m. heuristics_out[i] = Out-heuristic value of members[i]
  // Note: NOT NORMALIZED; MAY BE ALL-ZERO
  HeuristicsContainer heuristics_out;
  // List of size m. heuristics_out[i] = Out-heuristic value of members[i] as seed
  // Note: NOT NORMALIZED; MAY BE ALL-ZERO
  HeuristicsContainer heuristics_out_seed;
  // Matrix of shape (m, m). p_internal[i][j] = p of the directed edge members[i] -> members[j]
  PInternalContainer p_internal;
  // Matrix of shape (m, m). p_seed_internal[i][j] = p_seed of the directed edge members[i] -> members[j]
  PInternalContainer p_seed_internal;

  template <size_t N>
  auto get_members() const -> std::array<vertex_id_t, N> {
    BOOST_ASSERT_MSG(members.size() == N, "# of members mismatch.");
    auto res = std::array<vertex_id_t, N>{}; // Supports structured binding
    ranges::copy(members, res.begin());
    return res;
  }
};

struct CoarsenedEdgeDetails {
  static constexpr auto MAX_NUM_MEMBERS = CoarsenedVertexDetails::MAX_NUM_MEMBERS;
  using PCrossContainer = Array2D<edge_probability_t, MAX_NUM_MEMBERS, MAX_NUM_MEMBERS>;

  size_t n_members_left;
  size_t n_members_right;
  PCrossContainer p_cross;
  PCrossContainer p_seed_cross;
  WIMEdge merged;
};

struct CoarseningDetails {
  using VertexPair = std::pair<vertex_id_t, vertex_id_t>;
  using EdgeDetailsMap = std::map<VertexPair, CoarsenedEdgeDetails>;

  // N, # of vertices before coarsening
  vertex_id_t n;
  // Nc, # of vertices in the coarsened graph.
  vertex_id_t n_groups;
  // List of size N.
  std::vector<vertex_id_t> group_id;
  // List of size N.
  // Let {v1, v2, v3} be a group, with group_id[v1] = group_id[v2] = group_id[v3] = g, v1 < v2 < v3,
  // Then index_in_group[v1] = 0, index_in_group[v2] = 1, index_in_group[v3] = 2
  std::vector<vertex_id_t> index_in_group;
  // List of size Nc. groups[g] represents the group {v1, v2, v3} whose group index is g.
  std::vector<CoarsenedVertexDetails> groups;
  // [UNUSED] Map from vertex pair (u, v) to coarsened edge details, where 0 <= u, v < Nc
  // EdgeDetailsMap edges;

  template <std::convertible_to<vertex_id_t>... Args>
  auto to_group_id(Args... vertex_indices) const {
    BOOST_ASSERT_MSG(group_id.size() == n, "Size of group_id array mismatch with |V|.");
    if constexpr (sizeof...(Args) == 1) {
      return group_id[vertex_indices...];
    } else {
      return std::tuple{group_id[vertex_indices]...};
    }
  }

  template <std::convertible_to<vertex_id_t>... Args>
  auto to_group_ptr(Args... vertex_indices) const {
    BOOST_ASSERT_MSG(group_id.size() == n, "Size of group_id array mismatch with |V|.");
    if constexpr (sizeof...(Args) == 1) {
      return &groups[group_id[vertex_indices...]];
    } else {
      return std::tuple{&groups[group_id[vertex_indices]]...};
    }
  }

  template <std::convertible_to<vertex_id_t>... Args>
  auto to_index_in_group(Args... vertex_indices) const {
    return std::tuple{index_in_group[vertex_indices]...};
  }

  template <std::convertible_to<vertex_id_t>... Args>
  auto to_group_size(Args... vertex_indices) const {
    BOOST_ASSERT_MSG(group_id.size() == n, "Size of group_id array mismatch with |V|.");
    if constexpr (sizeof...(Args) == 1) {
      return groups[group_id[vertex_indices...]].members.size();
    } else {
      return std::tuple{groups[group_id[vertex_indices]].members.size()...};
    }
  }
};

struct CoarsenGraphResult {
  AdjacencyList<WIMEdge> coarsened_graph;
  InvAdjacencyList<WIMEdge> coarsened_inv_graph;
  std::vector<vertex_weight_t> coarsened_vertex_weights;
  CoarseningDetails details;

  auto in_degree(vertex_id_t v_coarsened) const {
    BOOST_ASSERT_MSG(v_coarsened >= 0 && v_coarsened < graph::num_vertices(coarsened_inv_graph),
                     "v is out of range [0, n) where n = # of coarsened vertices.");
    return graph::degree(coarsened_inv_graph, v_coarsened);
  }

  auto out_degree(vertex_id_t v_coarsened) const {
    BOOST_ASSERT_MSG(v_coarsened >= 0 && v_coarsened < graph::num_vertices(coarsened_graph),
                     "v is out of range [0, n) where n = # of coarsened vertices.");
    return graph::degree(coarsened_graph, v_coarsened);
  }

  auto degree(vertex_id_t v_coarsened) const {
    return in_degree(v_coarsened) + out_degree(v_coarsened);
  }
};

auto do_coarsen_wim_graph(const AdjacencyList<WIMEdge>& graph, const InvAdjacencyList<WIMEdge>& inv_graph,
                          vertex_id_t n_groups, std::span<const vertex_id_t> group_id,
                          const CoarseningParams& params) noexcept -> CoarsenGraphResult;

auto do_coarsen_wim_graph_w(const AdjacencyList<WIMEdge>& graph, const InvAdjacencyList<WIMEdge>& inv_graph,
                            std::span<const vertex_weight_t> vertex_weights, vertex_id_t n_groups,
                            std::span<const vertex_id_t> group_id, const CoarseningParams& params) noexcept
    -> CoarsenGraphResult;

// ---- Combination of Step 1 - 3 ----

auto coarsen_wim_graph(const AdjacencyList<WIMEdge>& graph, const InvAdjacencyList<WIMEdge>& inv_graph,
                       const CoarseningParams& params) noexcept -> CoarsenGraphResult;

auto coarsen_wim_graph_w(const AdjacencyList<WIMEdge>& graph, const InvAdjacencyList<WIMEdge>& inv_graph,
                         std::span<const vertex_weight_t> vertex_weights, const CoarseningParams& params) noexcept
    -> CoarsenGraphResult;

auto further_coarsen_wim_graph(const CoarsenGraphResult& last_result, const CoarseningParams& params) noexcept
    -> CoarsenGraphResult;

// ---- Step 4: Expanding seeds ----

struct SelectBestSeedResult {
  size_t index_in_group;
};

auto select_best_seed_in_group(const CoarsenedVertexDetails& v) noexcept -> SelectBestSeedResult;

struct ExpandSeedResult {
  std::vector<vertex_id_t> expanded_seeds;
};

auto expand_wim_seed_vertices(const AdjacencyList<WIMEdge>& graph, const CoarsenGraphResult& coarsening_result,
                              std::span<const vertex_id_t> coarsened_seeds, const ExpandingParams& params) noexcept
    -> rfl::Result<ExpandSeedResult>;

auto expand_wim_seed_vertices_w(const AdjacencyList<WIMEdge>& graph, std::span<const vertex_weight_t> vertex_weights,
                                const CoarsenGraphResult& coarsening_result,
                                std::span<const vertex_id_t> coarsened_seeds, const ExpandingParams& params) noexcept
    -> rfl::Result<ExpandSeedResult>;

inline auto further_expand_wim_seed_vertices(const CoarsenGraphResult& last_result,
                                             const CoarsenGraphResult& cur_result,
                                             std::span<const vertex_id_t> coarsened_seeds,
                                             const ExpandingParams& params) noexcept -> rfl::Result<ExpandSeedResult> {
  return expand_wim_seed_vertices_w(last_result.coarsened_graph, last_result.coarsened_vertex_weights, cur_result,
                                    coarsened_seeds, params);
}
