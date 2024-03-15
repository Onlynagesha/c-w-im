#pragma once

#include "utils/boost_container_common.h"
#include <boost/container/flat_set.hpp>

template <class T, class Comp = std::less<void>>
using FlatSet = boost::container::flat_set<T, Comp>;

template <class T, class Comp = std::less<void>>
using FlatMultiSet = boost::container::flat_multiset<T, Comp>;
