#pragma once

#include <rfl/internal/has_reflection_type_v.hpp>
#include <rfl/internal/is_rename.hpp>

namespace details {
template <class T, class = void>
struct ReflectionTypeHelper {
  using type = T;
};

template <class T>
struct ReflectionTypeHelper<T, std::enable_if_t<rfl::internal::has_reflection_type_v<T>>> {
  using type = typename T::ReflectionType;
};

template <class T>
struct ReflectionTypeHelper<T, std::enable_if_t<rfl::internal::is_rename_v<T>>> {
  using type = typename ReflectionTypeHelper<typename T::Type>::type;
};
} // namespace details

template <class T>
using ReflectionType = typename details::ReflectionTypeHelper<T>::type;

template <class Field>
using FieldReflectionType = ReflectionType<typename Field::Type>;

template <class T>
inline decltype(auto) reflected_value(const T& value) {
  if constexpr (rfl::internal::has_reflection_type_v<T>) {
    return value.value();
  } else {
    return value;
  }
}
