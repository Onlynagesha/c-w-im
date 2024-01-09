#pragma once

#include "create_dataset.h"
#include "graph_types.h"
#include "utils/utils.h"
#include "wim.h"
#include <fmt/format.h>

#define DUMP_REGISTERED_TYPES(F)          \
  F(WIMEdge)         /* graph_types.h */  \
  F(WBIMEdge)        /* graph_types.h */  \
  F(ReadGraphParams) /* read_dataset.h */ \
  F(WIMParams)       /* wim.h */

#define REGISTER_DUMP_FUNCTIONS(Type)                                   \
  auto dump(const Type& value, int indent = 0) noexcept -> std::string; \
  auto dump_as_json(const Type& value) noexcept -> json;

DUMP_REGISTERED_TYPES(REGISTER_DUMP_FUNCTIONS)

template <class T>
concept dump_methods_registered = requires(T t) {
  { dump(t) } -> std::same_as<std::string>;    // dump(T value)
  { dump(t, 4) } -> std::same_as<std::string>; // dump(T value, int indent)
};

template <dump_methods_registered T>
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
    auto dump_str = dump(value, indent);
    auto out = ctx.out();
    return ranges::copy(dump_str, out).out;
  }
};
