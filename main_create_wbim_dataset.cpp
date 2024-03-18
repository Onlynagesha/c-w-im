#include "main_create_dataset_common.h"

int main(int argc, char** argv) {
  return main_worker<WBIMEdge>(argc, argv)
      .and_then([](auto) {
        ELOG_INFO << "Done creating WBIM dataset.";
        return rfl::Result{0};
      })
      .or_else([](const rfl::Error& error) {
        ELOGFMT(CRITICAL, "ERROR: `{}'", error.what());
        return rfl::Result{-1};
      })
      .value_or(-1);
}
