#pragma once

#include <ylt/easylog.hpp>

#if defined(ENABLES_MYLOG_TRACE) && !defined(DISABLES_MYLOG_DEBUG)
  #define MYLOG_TRACE(...) ELOG_TRACE << (__VA_ARGS__)
  #define MYLOG_FMT_TRACE(...) ELOGFMT(TRACE, __VA_ARGS__)
#else
  #define MYLOG_TRACE(...) (void)0
  #define MYLOG_FMT_TRACE(...) (void)0
#endif

#if !defined(DISABLES_MYLOG_DEBUG)
  #define MYLOG_DEBUG(...) ELOG_DEBUG << (__VA_ARGS__)
  #define MYLOG_FMT_DEBUG(...) ELOGFMT(DEBUG, __VA_ARGS__)
#else
  #define MYLOG_DEBUG(...) (void)0
  #define MYLOG_FMT_DEBUG(...) (void)0
#endif
