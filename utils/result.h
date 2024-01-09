#pragma once

#include "utils/demangle.h"
#include <rfl/Result.hpp>
#include <stacktrace>

constexpr auto RESULT_VOID_SUCCESS = rfl::Nothing{};

using ResultVoid = rfl::Result<rfl::Nothing>;

#define RFL_RESULT_CATCH_HANDLER()                                                                               \
  catch (std::exception & e) {                                                                                   \
    constexpr auto msg_pattern = "Exception of type {} caught: `{}'\nStacktrace where exception is caught:\n{}"; \
    auto st = to_string(std::stacktrace::current());                                                             \
    return rfl::Error{fmt::format(msg_pattern, demangle_type_name(e), e.what(), st)};                            \
  }                                                                                                              \
  catch (...) {                                                                                                  \
    constexpr auto msg_pattern = "Unknown exception caught. Stacktrace where exception is caught:\n{}";          \
    auto st = to_string(std::stacktrace::current());                                                             \
    return rfl::Error{fmt::format(msg_pattern, st)};                                                             \
  }
