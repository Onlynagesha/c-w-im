#pragma once

#include <nwgraph/util/demangle.hpp>

template <class T>
inline auto demangle_type_name() {
  return nw::graph::demangle(typeid(T).name(), nullptr, nullptr, nullptr);
}

template <class T>
inline auto demangle_type_name(const T&) {
  return demangle_type_name<T>();
}
