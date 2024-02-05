#include <chrono>
#include <expected>
#include <fmt/color.h>
#include <iostream>
#include <rfl/Result.hpp>

struct VeryLarge {
  static inline auto default_ctor_count = 0;

  VeryLarge() : value(std::to_string(default_ctor_count++)) {
    std::cout << "Default constructor: VeryLarge()\n";
  }
  VeryLarge(VeryLarge&& other) noexcept : value(std::move(other.value)) {
    std::cout << fmt::format(fg(fmt::color::yellow), "Move constructor: VeryLarge(VeryLarge&&)\n");
    other.value = "(Moved away)";
  }
  VeryLarge(const VeryLarge& other) : value(other.value) {
    std::cout << fmt::format(fg(fmt::color::orange), "Copy constructor: VeryLarge(const VeryLarge&)\n");
  }

  std::string value;
};

template <class T>
using ResultType = rfl::Result<T>;
// using ResultType = std::expected<T, rfl::Error>;

auto some_func() -> ResultType<VeryLarge> {
  return VeryLarge{};
}

int main() {
  auto foo = some_func();
  std::cout << fmt::format("Before foo.and_then(): foo.value = {:?}\n", (*foo).value);
  // Here we attempt to modify the contents inside foo via a left-reference.
  foo.and_then([](VeryLarge& obj) -> ResultType<rfl::Nothing> {
    obj.value += "-modified"; // Modification here
    std::cout << fmt::format("foo.and_then(): obj.value: {:?}\n", obj.value);
    return rfl::Nothing{};
  });
  // We expect that the modification above is represented here.
  std::cout << fmt::format("After foo.and_then(): foo.value: {:?}\n", (*foo).value);

  // Still, for right-values, copy constructor shall not be triggered.
  some_func().and_then([](VeryLarge obj) -> ResultType<rfl::Nothing> {
    std::cout << fmt::format("some_func().and_then(): obj.value: {:?}\n", obj.value);
    return rfl::Nothing{};
  });

  return 0;
}
