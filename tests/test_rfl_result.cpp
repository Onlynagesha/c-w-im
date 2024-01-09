#include <chrono>
#include <fmt/color.h>
#include <fmt/format.h>
#include <iostream>
#include <rfl/Result.hpp>
#include <stacktrace>
#include <thread>

void simulate_time_consuming(int seconds) {
  constexpr auto msg_pattern = "Simulating very time-consuming operation: sleeps for {} sec.\n";
  std::cout << fmt::format(fg(fmt::color::orange_red), msg_pattern, seconds);
  if (seconds > 0) {
    std::this_thread::sleep_for(std::chrono::seconds{seconds});
  }
}

struct VeryLarge {
  static inline int default_ctor_count = 0;
  static inline int copy_ctor_count = 0;
  static inline int move_ctor_count = 0;

  VeryLarge() {
    default_ctor_count += 1;

    auto st = to_string(std::stacktrace::current(1));
    std::cout << fmt::format("Default constructor: VeryLarge()\nStacktrace:\n{}", st);
  }

  VeryLarge(const VeryLarge&) {
    copy_ctor_count += 1;

    constexpr auto msg_pattern = "Copy constructor: VeryLarge(const VeryLarge&)\nStacktrace:\n{}";
    auto st = to_string(std::stacktrace::current(1));
    std::cout << fmt::format(fg(fmt::color::yellow), msg_pattern, st);
    simulate_time_consuming(3);
  }

  VeryLarge(VeryLarge&& src) noexcept {
    move_ctor_count += 1;

    constexpr auto msg_pattern = "Move constructor: VeryLarge(VeryLarge&&)\nStacktrace:\n{}";
    auto st = to_string(std::stacktrace::current(1));
    std::cout << fmt::format(fg(fmt::color::light_green), msg_pattern, st);
  }
};

int main() {
  auto result1 = rfl::Result<VeryLarge>{VeryLarge{}};
  result1.and_then([](VeryLarge&) -> rfl::Result<int> {
    std::cout << fmt::format("result1.and_then(), with result1 as non-const left-reference.\n");
    return 0;
  });

  const auto result2 = rfl::Result<VeryLarge>{VeryLarge{}};
  result2.and_then([](const VeryLarge&) -> rfl::Result<int> {
    std::cout << fmt::format("result2.and_then(), with result2 as const left-reference.\n");
    return 0;
  });

  std::move(result1).and_then([](VeryLarge) -> rfl::Result<int> {
    std::cout << fmt::format("move(result1).and_then(), with result1 as non-const right-reference.\n");
    return 0;
  });

  std::move(result2).and_then([](const VeryLarge&) -> rfl::Result<int> {
    std::cout << fmt::format("move(result2).and_then(), with result2 as const right-reference.\n");
    return 0;
  });

  return 0;
}
