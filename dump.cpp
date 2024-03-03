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
} // namespace

#define IMPLEMENT_DUMP_FUNCTIONS_AUTO(Type)                          \
  auto dump(const Type& value, int indent) noexcept -> std::string { \
    return ::dump_generic(value, indent);                            \
  }                                                                  \
  auto dump_as_json(const Type& value) noexcept -> json {            \
    return ::dump_as_json_generic(value);                            \
  }

DUMP_REGISTERED_TYPES_AUTO_GENERATED(IMPLEMENT_DUMP_FUNCTIONS_AUTO)
#undef IMPLEMENT_DUMP_FUNCTIONS_AUTO
