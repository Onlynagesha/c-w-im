#pragma once

// HINT: It's recommended to include this "heavy" header
// (which involves huge amount of metaprogramming operations)
// into some specific source file instead of using it as header-only
// to reduce compilation overhead.

#include "utils/reflect.h"
#include "utils/utils.h"
#include <nlohmann/json.hpp>
#include <rfl/named_tuple_t.hpp>
#include <rfl/to_view.hpp>

namespace details {
template <class Field, class TView, class Func>
inline auto for_each_field_visitor(const TView& value_view, Func&& func) {
  constexpr auto name_literal = Field::name_;
  std::invoke(func, name_literal.str(), *value_view.template get<name_literal>());
}

template <class FieldTuple, class TView, class Func, size_t... Indices>
inline auto for_each_field_impl(const TView& value_view, Func&& func, std::index_sequence<Indices...>) -> void {
  (for_each_field_visitor<std::tuple_element_t<Indices, FieldTuple>>(value_view, func), ...);
}

template <class FieldTuple, class T, class Func>
inline auto for_each_field(T& value, Func&& func) -> void {
  auto value_view = rfl::to_view(value);
  for_each_field_impl<FieldTuple>(value_view, func, std::make_index_sequence<std::tuple_size_v<FieldTuple>>{});
}
} // namespace details

template <class T, class Func>
inline auto for_each_field(T& value, Func&& func) -> void {
  using FieldTuple = typename rfl::named_tuple_t<T>::Fields;
  details::for_each_field<FieldTuple>(value, func);
}

template <class T>
inline auto read_fields_from_json(T& value, const json& json_root) {
  for_each_field(value, [&json_root]<class Field>(std::string_view key, Field& field) {
    auto it = json_root.find(key);
    if (it != json_root.end()) {
      // Triggers check during operator= for rfl::Validator types
      field = it->template get<ReflectionType<Field>>();
    }
  });
}
