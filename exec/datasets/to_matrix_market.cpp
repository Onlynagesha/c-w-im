#include "to_matrix_market.h"
#include <ylt/easylog.hpp>

auto main_worker(int argc, char** argv) -> ResultVoid try {
  return ToMatrixMarketParams::parse_from_args(argc, argv).and_then([](ToMatrixMarketParams params) {
    return graph_text_to_matrix_market(params);
  });
}
RFL_RESULT_CATCH_HANDLER()

int main(int argc, char** argv) {
  easylog::set_min_severity(easylog::Severity::DEBUG);
  return main_worker(argc, argv)
      .transform([](auto) {
        ELOG_INFO << "Done transforming graph text to Matrix Market format.";
        return 0;
      })
      .or_else([](const rfl::Error& error) {
        ELOGFMT(CRITICAL, "Failed transforming graph text to Matrix Market format: `{}'", error.what());
        return -1;
      })
      .value_or(-1);
}
