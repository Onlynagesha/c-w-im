#pragma once

#include <boost/container/static_vector.hpp>

template <class T, size_t Capacity>
using StaticVector = boost::container::static_vector<T, Capacity>;