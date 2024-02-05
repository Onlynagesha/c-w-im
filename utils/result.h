#pragma once

#include "utils/demangle.h"
#include <boost/stacktrace.hpp>
#include <execinfo.h>
#include <fmt/format.h>
#include <rfl/Result.hpp>

constexpr auto RESULT_VOID_SUCCESS = rfl::Nothing{};

using ResultVoid = rfl::Result<rfl::Nothing>;

#define RFL_RESULT_CATCH_HANDLER()                                                                    \
  catch (std::exception & e) {                                                                        \
    constexpr auto msg_pattern = "Exception of type {} caught at line #{} of file {}: `{}'";          \
    return rfl::Error{fmt::format(msg_pattern, demangle_type_name(e), __LINE__, __FILE__, e.what())}; \
  }                                                                                                   \
  catch (...) {                                                                                       \
    constexpr auto msg_pattern = "Unknown exception caught at line #{} of file {}.";                  \
    return rfl::Error{fmt::format(msg_pattern, __LINE__, __FILE__)};                                  \
  }
