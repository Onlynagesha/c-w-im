#pragma once

#include "coarsening.h"
#include "create_dataset.h"
#include "experiments.h"
#include "graph_types.h"
#include "utils/utils.h"
#include <fmt/format.h>

#define DUMP_REGISTERED_TYPES_AUTO_GENERATED(F)          \
  F(WIMEdge)                        /* graph_types.h */  \
  F(WBIMEdge)                       /* graph_types.h */  \
  F(CreateWIMDatasetParams)         /* read_dataset.h */ \
  F(CreateWBIMDatasetParams)        /* read_dataset.h */ \
  F(CoarseningParams)               /* coarsening.h */   \
  F(WIMSketchingExperimentParams)   /* experiments.h */  \
  F(WBIMSketchingExperimentParams)  /* experiments.h */  \
  F(WIMCoarseningExperimentParams)  /* experiments.h */  \
  F(WBIMCoarseningExperimentParams) /* experiments.h */  \
  F(WIMContrastExperimentParams)    /* experiments.h */

#define REGISTER_AUTO_GENERATED_DUMP_FUNCTIONS_WITH_JSON(Type)          \
  auto dump(const Type& value, int indent = 0) noexcept -> std::string; \
  auto dump_as_json(const Type& value) noexcept -> json;

DUMP_REGISTERED_TYPES_AUTO_GENERATED(REGISTER_AUTO_GENERATED_DUMP_FUNCTIONS_WITH_JSON)
#undef REGISTER_AUTO_GENERATED_DUMP_FUNCTIONS_WITH_JSON

template <class T>
concept has_dump_free_methods = requires(const T& t) {
  { dump(t) } -> std::same_as<std::string>;
  { dump(t, 4) } -> std::same_as<std::string>;
};

template <class T>
concept has_dump_member_methods = requires(const T& t) {
  { t.dump() } -> std::same_as<std::string>;
  { t.dump(4) } -> std::same_as<std::string>;
};

template <class T>
concept has_dump_methods = has_dump_free_methods<T> || has_dump_member_methods<T>;

template <has_dump_methods T>
struct fmt::formatter<T> {
  int indent = 0;

  constexpr auto parse(fmt::format_parse_context& ctx) {
    auto it = ctx.begin();
    for (; it != ctx.end() && *it != '}'; ++it) {
      if (*it >= '0' && *it <= '9') {
        indent = indent * 10 + (*it - '0');
      } else {
        throw_format_error("Invalid format specifier: An decimal non-negative integer expected.");
      }
    }
    return it;
  }

  auto format(const T& value, fmt::format_context& ctx) const {
    auto dump_str = [&]() {
      if constexpr (has_dump_free_methods<T>) {
        return dump(value, indent);
      } else {
        return value.dump(indent);
      }
    }();
    auto out = ctx.out();
    return ranges::copy(dump_str, out).out;
  }
};
