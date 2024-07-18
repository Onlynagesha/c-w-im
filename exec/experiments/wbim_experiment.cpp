#include "experiments.h"

int main(int argc, char** argv) {
  return wbim_experiment(argc, argv)
      .and_then([](auto) {
        ELOG_INFO << "WBIM experiment done.";
        return rfl::Result{0};
      })
      .or_else([](const rfl::Error& error) {
        ELOGFMT(CRITICAL, "WBIM experiment halts due to error: `{}'", error.what());
        return rfl::Result{-1};
      })
      .value_or(-1);
}
