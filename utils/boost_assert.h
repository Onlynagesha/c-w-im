#pragma once

#include <boost/assert.hpp>
#include <execinfo.h>
#include <fmt/format.h>

namespace details {
inline void dump_assert_fail(char const* expr, char const* msg, char const* function, char const* file, long line) {
  constexpr auto msg_pattern = "Assertion failed: {0}"
                               "\n\tDescription: {1}"
                               "\n\tSource location: {2}:{3}"
                               "\n\tAt function: {4}";
  fmt::print(stderr, msg_pattern, expr, msg, file, line, function);
}
} // namespace details

namespace boost {
inline void assertion_failed(char const* expr, char const* function, char const* file, long line) {
  details::dump_assert_fail(expr, "(None)", function, file, line);
  std::abort();
}
inline void assertion_failed_msg(char const* expr, char const* msg, char const* function, char const* file, long line) {
  details::dump_assert_fail(expr, msg, function, file, line);
  std::abort();
}
} // namespace boost
