#pragma once

#include "graph_types.h"
#include "utils/static_vector.h"
#include "vertex_set.h"

#define COARSENING_DETAILS_TYPES(F) \
  F(WIMCoarsenedVertexBrief)        \
  F(WIMCoarsenedVertexDetails)      \
  F(WBIMCoarsenedVertexBrief)       \
  F(WBIMCoarsenedVertexDetails)     \
  F(WIMCoarsenedEdgeDetails)        \
  F(WBIMCoarsenedEdgeDetails)       \
  F(WIMCoarseningBrief)             \
  F(WIMCoarseningDetails)           \
  F(WBIMCoarseningBrief)            \
  F(WBIMCoarseningDetails)          \
  F(WIMCoarsenGraphBriefResult)     \
  F(WIMCoarsenGraphDetailedResult)  \
  F(WBIMCoarsenGraphBriefResult)    \
  F(WBIMCoarsenGraphDetailedResult)

// ---- COARSENED VERTEX DETAILS ----

struct CoarsenedVertexBriefBase {
  // Mongoose algorithm ensures that the size of each group is at most 3 (match + 1 adopted)
  static constexpr auto MAX_N_MEMBERS = size_t{3};
  using MemberContainer = StaticVector<vertex_id_t, MAX_N_MEMBERS>;

  // List of size M, # of members in the groups, i.e. the vertices before coarsening into current group
  MemberContainer members = {};

  CoarsenedVertexBriefBase() = default;
  explicit CoarsenedVertexBriefBase(MemberContainer members) : members(std::move(members)) {}

  auto n_members() const -> size_t {
    return members.size();
  }

  // Equivalent to range(n_members()), used to reduce one layer of parentheses :)
  auto member_indices() const -> ranges::iota_view<size_t, size_t> {
    return range(n_members());
  }
};

struct WIMCoarsenedVertexBrief : public CoarsenedVertexBriefBase {
  // Index in range [0, M), the best expanded seed by local information.
  size_t best_seed_index;

  WIMCoarsenedVertexBrief() = default;
  WIMCoarsenedVertexBrief(MemberContainer members, size_t best_seed_index)
      : CoarsenedVertexBriefBase(std::move(members)), best_seed_index(best_seed_index) {}

  auto dump(int indent = 0, int level = 0) const noexcept -> std::string;
};

struct WBIMCoarsenedVertexBrief : public CoarsenedVertexBriefBase {
  // Index in range [0, M), the best expanded boosted vertex by local information
  size_t best_boosted_index;

  WBIMCoarsenedVertexBrief() = default;
  WBIMCoarsenedVertexBrief(MemberContainer members, size_t best_boosted_index)
      : CoarsenedVertexBriefBase(std::move(members)), best_boosted_index(best_boosted_index) {}

  auto dump(int indent = 0, int level = 0) const noexcept -> std::string;
};

struct CoarsenedVertexDetailsBase : public CoarsenedVertexBriefBase {
  using VertexWeightsContainer = StaticVector<vertex_weight_t, MAX_N_MEMBERS>;
  using HeuristicsContainer = StaticVector<edge_probability_t, MAX_N_MEMBERS>;
  using PInternalContainer = Array2D<edge_probability_t, MAX_N_MEMBERS, MAX_N_MEMBERS>;

  // List of size M. vertex_weights[i] = Vertex weight of members[i]
  VertexWeightsContainer vertex_weights = {};
  // List of size M. heuristics_in[i] = In-heuristic value of members[i]
  // Note: In range [0, +inf); MAY BE ALL-ZERO
  HeuristicsContainer heuristics_in = {};
  // List of size M. heuristics_out[i] = Out-heuristic value of members[i]
  // Note: In range [0, +inf); MAY BE ALL-ZERO
  HeuristicsContainer heuristics_out = {};
  // Matrix of shape (M, M). p_internal[i][j] = p of the directed edge members[i] -> members[j]
  PInternalContainer p_internal = {};
};

struct WIMCoarsenedVertexDetails : public CoarsenedVertexDetailsBase {
  using BriefType = WIMCoarsenedVertexBrief;

  // List of size M. heuristics_out[i] = Out-heuristic value of members[i] as seed
  // Note: In range [0, +inf); MAY BE ALL-ZERO
  HeuristicsContainer heuristics_out_seed = {};
  // Matrix of shape (M, M). p_seed_internal[i][j] = p_seed of the directed edge members[i] -> members[j]
  PInternalContainer p_seed_internal = {};
  // Index in range [0, M), the best expanded seed by local information.
  size_t best_seed_index;

  auto to_brief() const -> WIMCoarsenedVertexBrief {
    return {members, best_seed_index};
  }

  auto dump(int indent = 0, int level = 0) const noexcept -> std::string;
};

struct WBIMCoarsenedVertexDetails : public CoarsenedVertexDetailsBase {
  using BriefType = WBIMCoarsenedVertexBrief;

  bool is_seed = false;
  // List of size M. heuristics_in_boost[i] = In-heuristic value of members[i] assuming it's boosted
  // Note: In range [0, +inf); MAY BE ALL-ZERO
  HeuristicsContainer heuristics_in_boost = {};
  // Matrix of shape (M, M). p_boost_internal[i][j] = p_boost of the directed edge members[i] -> members[j]
  PInternalContainer p_boost_internal = {};
  // Index in range [0, M), the best expanded boosted vertex by local information
  size_t best_boosted_index;

  auto to_brief() const -> WBIMCoarsenedVertexBrief {
    return {members, best_boosted_index};
  }

  auto dump(int indent = 0, int level = 0) const noexcept -> std::string;
};

template <class T>
concept is_coarsened_vertex_brief = same_as_either<T, WIMCoarsenedVertexBrief, WBIMCoarsenedVertexBrief>;
template <class T>
concept is_coarsened_vertex_details = same_as_either<T, WIMCoarsenedVertexDetails, WBIMCoarsenedVertexDetails>;
template <class T>
concept is_wim_coarsened_vertex_info = same_as_either<T, WIMCoarsenedVertexBrief, WBIMCoarsenedVertexBrief>;
template <class T>
concept is_wbim_coarsened_vertex_info = same_as_either<T, WIMCoarsenedVertexDetails, WBIMCoarsenedVertexDetails>;
template <class T>
concept is_coarsened_vertex_info = is_wim_coarsened_vertex_info<T> || is_wbim_coarsened_vertex_info<T>;

// ---- COARSENED EDGE DETAILS ----

/*
Details of the coarsened edge Gu -> Gv,
where the coarsened vertex Gu -> {u1, u2, ...}, Gv -> {v1, v2, ...} after expansion.
*/
struct CoarsenedEdgeDetailsBase {
  static constexpr auto MAX_N_MEMBERS = WIMCoarsenedVertexDetails::MAX_N_MEMBERS;
  using PCrossContainer = Array2D<edge_probability_t, MAX_N_MEMBERS, MAX_N_MEMBERS>;

  // Size of {u1, u2, ...}. Demoted as Mu
  size_t n_members_left;
  // Size of {v1, v2, ...}. Denoted as Mv
  size_t n_members_right;
  // p_cross[i][j] = p(ui, vj)
  PCrossContainer p_cross;

  CoarsenedEdgeDetailsBase() = default;
  CoarsenedEdgeDetailsBase(size_t n_members_left, size_t n_members_right)
      : n_members_left(n_members_left), n_members_right(n_members_right), p_cross() {}
};

struct WIMCoarsenedEdgeDetails : public CoarsenedEdgeDetailsBase {
  // p_seed_cross[i][j] = p_seed(ui, vj)
  PCrossContainer p_seed_cross;
  // Properties of the merged edge
  WIMEdge merged;

  WIMCoarsenedEdgeDetails() = default;
  WIMCoarsenedEdgeDetails(size_t n_members_left, size_t n_members_right)
      : CoarsenedEdgeDetailsBase(n_members_left, n_members_right), p_seed_cross(), merged() {}

  auto dump(int indent = 0, int level = 0) const noexcept -> std::string;
};

struct WBIMCoarsenedEdgeDetails : public CoarsenedEdgeDetailsBase {
  // p_boost_cross[i][j] = p_boost(ui, vj)
  PCrossContainer p_boost_cross;
  // Properties of the merged edge
  WBIMEdge merged;

  WBIMCoarsenedEdgeDetails() = default;
  WBIMCoarsenedEdgeDetails(size_t n_members_left, size_t n_members_right)
      : CoarsenedEdgeDetailsBase(n_members_left, n_members_right), p_boost_cross(), merged() {}

  auto dump(int indent = 0, int level = 0) const noexcept -> std::string;
};

template <class T>
concept is_coarsened_edge_details = same_as_either<T, WIMCoarsenedEdgeDetails, WBIMCoarsenedEdgeDetails>;

// ---- COARSENING INFO (OF COARSENED GROUPS) ----

struct CoarseningInfoTags {
  struct value_initialize_groups_t {};
  static constexpr auto value_initialize_groups = value_initialize_groups_t{};

  struct group_as_empty_list_t {};
  static constexpr auto group_as_empty_list = group_as_empty_list_t{};
};

template <is_coarsened_vertex_info VertexInfoType>
struct CoarseningInfoBase : public CoarseningInfoTags {
  // N, # of vertices before coarsening
  vertex_id_t n;
  // Nc, # of vertices in the coarsened graph.
  vertex_id_t n_coarsened;
  // List of size Nc. groups[g] represents the group {v1, v2, ...} whose group index is g,
  // i.e. all the vertices v1, v2, ... that are coarsened to g.
  std::vector<VertexInfoType> groups;

  CoarseningInfoBase(group_as_empty_list_t, vertex_id_t n, vertex_id_t n_coarsened) : n(n), n_coarsened(n_coarsened) {}

  // Creates a value-initialized group list of size Nc.
  CoarseningInfoBase(value_initialize_groups_t, vertex_id_t n, vertex_id_t n_coarsened)
      : n(n), n_coarsened(n_coarsened), groups(n_coarsened) {}
};

template <is_coarsened_vertex_brief VertexBriefType>
struct CoarseningBrief : public CoarseningInfoBase<VertexBriefType> {
  // Inherited constructors
  using CoarseningInfoBase<VertexBriefType>::CoarseningInfoBase;

  auto dump(int indent = 0, int level = 0) const noexcept -> std::string;
};

using WIMCoarseningBrief = CoarseningBrief<WIMCoarsenedVertexBrief>;
using WBIMCoarseningBrief = CoarseningBrief<WBIMCoarsenedVertexBrief>;

template <is_coarsened_vertex_details VertexDetailsType>
struct CoarseningDetails : public CoarseningInfoBase<VertexDetailsType> {
  using BaseType = CoarseningInfoBase<VertexDetailsType>;
  using VertexBriefType = typename VertexDetailsType::BriefType;

  using group_as_empty_list_t = CoarseningInfoTags::group_as_empty_list_t;
  using value_initialize_groups_t = CoarseningInfoTags::value_initialize_groups_t;

  using BaseType::groups;

  // List of size N.
  std::vector<vertex_id_t> group_id;
  // List of size N.
  // Let {v1, v2, v3} be a group, with group_id[v1] = group_id[v2] = group_id[v3] = g, v1 < v2 < v3,
  // Then index_in_group[v1] = 0, index_in_group[v2] = 1, index_in_group[v3] = 2
  std::vector<vertex_id_t> index_in_group;

  template <same_as_either<group_as_empty_list_t, value_initialize_groups_t> Tag>
  CoarseningDetails(Tag tag, vertex_id_t n, vertex_id_t n_coarsened, std::vector<vertex_id_t> group_id,
                    std::vector<vertex_id_t> index_in_group)
      : BaseType(tag, n, n_coarsened), group_id(std::move(group_id)), index_in_group(std::move(index_in_group)) {}

  auto to_brief() const -> CoarseningBrief<VertexBriefType> {
    auto res = CoarseningBrief<VertexBriefType>(group_as_empty_list_t{}, this->n, this->n_coarsened);
    for (const auto& g : this->groups) {
      res.groups.push_back(g.to_brief());
    }
    return res;
  }

  auto dump(int indent = 0, int level = 0) const noexcept -> std::string;

#define COARSENING_DETAILS_ACCESSOR_FUNCTIONS(F)                     \
  F(to_group_id, vertex_indices, group_id[_1])                       \
  F(to_group_ptr, vertex_indices, groups.data() + group_id[_1])      \
  F(to_index_in_group, vertex_indices, index_in_group[_1])           \
  F(to_group_size, vertex_indices, groups[group_id[_1]].n_members()) \
  F(group_id_to_size, group_indices, groups[_1].n_members())

#define DEFINE_COARSENING_DETAILS_ACCESSOR_FUNCTION(function_name, args_name, expr) \
  template <std::convertible_to<vertex_id_t>... Args>                               \
  auto function_name(Args... args_name) const {                                     \
    auto get_fn = LAMBDA_1(expr);                                                   \
    if constexpr (sizeof...(Args) == 1) {                                           \
      return get_fn(args_name...);                                                  \
    } else {                                                                        \
      return std::tuple{get_fn(args_name)...};                                      \
    }                                                                               \
  }
  // Accessors
  COARSENING_DETAILS_ACCESSOR_FUNCTIONS(DEFINE_COARSENING_DETAILS_ACCESSOR_FUNCTION)

#undef COARSENING_DETAILS_ACCESSOR_FUNCTIONS
#undef DEFINE_COARSENING_DETAILS_ACCESSOR_FUNCTION
};

using WIMCoarseningDetails = CoarseningDetails<WIMCoarsenedVertexDetails>;
using WBIMCoarseningDetails = CoarseningDetails<WBIMCoarsenedVertexDetails>;

// ---- COARSENING RESULT ----

template <is_edge_property E, class BriefOrDetailsType>
struct CoarsenGraphResultBase {
  AdjacencyListPair<E> coarsened;
  BriefOrDetailsType details;

  CoarsenGraphResultBase(AdjacencyListPair<E> coarsened, BriefOrDetailsType details)
      : coarsened(std::move(coarsened)), details(std::move(details)) {}

  CoarsenGraphResultBase(const DirectedEdgeList<E>& coarsened_edges, //
                         std::vector<vertex_weight_t> vertex_weights, BriefOrDetailsType details)
      : details(std::move(details)) {
    coarsened = {
        .adj_list = {as_non_const(coarsened_edges)},
        .inv_adj_list = {as_non_const(coarsened_edges)},
        .vertex_weights = std::move(vertex_weights),
    };
  }
};

template <class BriefOrDetailsType>
struct WIMCoarsenGraphResult : public CoarsenGraphResultBase<WIMEdge, BriefOrDetailsType> {
  // Inherited constructors
  using CoarsenGraphResultBase<WIMEdge, BriefOrDetailsType>::CoarsenGraphResultBase;

  auto dump(int indent = 0, int level = 0) const noexcept -> std::string;
};

using WIMCoarsenGraphBriefResult = WIMCoarsenGraphResult<WIMCoarseningBrief>;
using WIMCoarsenGraphDetailedResult = WIMCoarsenGraphResult<WIMCoarseningDetails>;

template <class BriefOrDetailsType>
struct WBIMCoarsenGraphResult : public CoarsenGraphResultBase<WBIMEdge, BriefOrDetailsType> {
  using BaseType = CoarsenGraphResultBase<WBIMEdge, BriefOrDetailsType>;
  // Maps seeds to coarsened index
  VertexSet coarsened_seeds;

  WBIMCoarsenGraphResult(WBIMAdjacencyListPair coarsened, VertexSet coarsened_seeds, BriefOrDetailsType details)
      : BaseType(std::move(coarsened), std::move(details)), coarsened_seeds(std::move(coarsened_seeds)) {}

  WBIMCoarsenGraphResult(const DirectedEdgeList<WBIMEdge>& coarsened_edges, std::vector<vertex_weight_t> vertex_weights,
                         VertexSet coarsened_seeds, BriefOrDetailsType details)
      : BaseType(coarsened_edges, std::move(vertex_weights), std::move(details)),
        coarsened_seeds(std::move(coarsened_seeds)) {}

  auto dump(int indent = 0, int level = 0) const noexcept -> std::string;
};

using WBIMCoarsenGraphBriefResult = WBIMCoarsenGraphResult<WBIMCoarseningBrief>;
using WBIMCoarsenGraphDetailedResult = WBIMCoarsenGraphResult<WBIMCoarseningDetails>;

// ---- TYPE TRAITS ----

template <is_edge_property E>
struct CoarseningDataTraits {
  static_assert(rfl::always_false_v<E>, "Traits does not exist.");
};

template <>
struct CoarseningDataTraits<WIMEdge> {
  using VertexBrief = WIMCoarsenedVertexBrief;
  using VertexDetails = WIMCoarsenedVertexDetails;
  using EdgeDetails = WIMCoarsenedEdgeDetails;
  using CoarseningBrief = WIMCoarseningBrief;
  using CoarseningDetails = WIMCoarseningDetails;
  using BriefResult = WIMCoarsenGraphBriefResult;
  using DetailedResult = WIMCoarsenGraphDetailedResult;
};

template <>
struct CoarseningDataTraits<WBIMEdge> {
  using VertexBrief = WBIMCoarsenedVertexBrief;
  using VertexDetails = WBIMCoarsenedVertexDetails;
  using EdgeDetails = WBIMCoarsenedEdgeDetails;
  using CoarseningBrief = WBIMCoarseningBrief;
  using CoarseningDetails = WBIMCoarseningDetails;
  using BriefResult = WBIMCoarsenGraphBriefResult;
  using DetailedResult = WBIMCoarsenGraphDetailedResult;
};

// ---- DUMP FUNCTIONS ----

// Note: declaration is required to trigger compilation in coarsening_dump.cpp
#define DECLARE_FREE_DUMP_FUNCTION_FOR_COARSENING_DETAILS(Type) \
  auto dump(const Type& obj, int indent = 0, int level = 0) noexcept -> std::string;

COARSENING_DETAILS_TYPES(DECLARE_FREE_DUMP_FUNCTION_FOR_COARSENING_DETAILS)
#undef DECLARE_FREE_DUMP_FUNCTION_FOR_COARSENING_DETAILS
