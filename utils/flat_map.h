#pragma once

#include "utils/boost_container_common.h"
#include <boost/container/flat_map.hpp>

template <class Key, class Value, class Comp = std::less<void>>
using FlatMap = boost::container::flat_map<Key, Value, Comp>;

template <class Key, class Value, class Comp = std::less<void>>
using FlatMultiMap = boost::container::flat_multimap<Key, Value, Comp>;
