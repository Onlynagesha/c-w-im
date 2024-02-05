#pragma once

#include <concepts>
#include <functional>
#include <iterator>
#include <nlohmann/json_fwd.hpp>
#include <numeric>
#include <random>
#include <ranges>
#include <utility>

#define LAMBDA_0(...) [&]() { return (__VA_ARGS__); }
#define LAMBDA_1(...) [&](auto&& _1) { return (__VA_ARGS__); }
#define LAMBDA_2(...) [&](auto&& _1, auto&& _2) { return (__VA_ARGS__); }

#ifdef NDEBUG
  #define NOEXCEPT_IF_NDEBUG noexcept
#else
  #define NOEXCEPT_IF_NDEBUG
#endif

using namespace std::literals;

namespace ranges = std::ranges;
namespace views = std::views;

using nlohmann::json;

using ssize_t = std::make_signed_t<size_t>;

// Alternative of std::forward_like: https://en.cppreference.com/w/cpp/utility/forward_like
template <class T, class U>
[[nodiscard]] constexpr auto&& forward_like(U&& x) noexcept {
  constexpr bool is_adding_const = std::is_const_v<std::remove_reference_t<T>>;
  if constexpr (std::is_lvalue_reference_v<T&&>) {
    if constexpr (is_adding_const) {
      return std::as_const(x);
    } else {
      return static_cast<U&>(x);
    }
  } else {
    if constexpr (is_adding_const) {
      return std::move(std::as_const(x));
    } else {
      return std::move(x);
    }
  }
}

template <class T>
constexpr auto is_std_vector_v = false;

template <class T, class Alloc>
constexpr auto is_std_vector_v<std::vector<T, Alloc>> = true;

template <class T, class... Candidates>
concept same_as_either = (std::same_as<T, Candidates> || ...);

template <class T, class... Candidates>
concept convertible_to_either = (std::convertible_to<T, Candidates> || ...);

template <class Func, class... Args>
concept invocable_or_nullptr = std::same_as<Func, std::nullptr_t> || std::invocable<Func, Args...>;

template <class Func, class Ret, class... Args>
concept invocable_r_or_nullptr = std::same_as<Func, std::nullptr_t> || std::is_invocable_r_v<Ret, Func, Args...>;

// WARNING: No thread safety.
inline auto rand_engine = std::minstd_rand{std::random_device{}()};

inline auto rand_float() -> float {
  constexpr auto ieee754_float_fraction_bits = 23u;
  constexpr auto ieee754_float_fraction_mask = (1u << ieee754_float_fraction_bits) - 1;
  constexpr auto ieee754_float_non_fraction_part = 127u << ieee754_float_fraction_bits;

  auto u = rand_engine();
  // res as float32 = 1.??? x 2^0 = random value in range [1, 2)
  auto res = (u & ieee754_float_fraction_mask) | ieee754_float_non_fraction_part;
  return reinterpret_cast<float&>(res) - 1;
}

inline auto rand_bool(float p) -> bool {
  return rand_float() < p;
}

template <std::integral IntType>
inline auto range(IntType n) {
  return views::iota(static_cast<IntType>(0), n);
}

template <std::integral IntType>
inline auto range(IntType first, IntType last) {
  return views::iota(first, last);
}

template <std::input_or_output_iterator Iter, std::sentinel_for<Iter> Sentinel>
inline auto range(Iter first, Sentinel last) {
  return ranges::subrange{first, last};
}

template <class T, ranges::input_range Range>
  requires(std::is_convertible_v<ranges::range_value_t<Range>, T>)
inline auto accumulate_sum(Range&& range, T init = T{}) -> T {
  for (auto&& elem : range) {
    init = init + forward_like<Range>(elem);
  }
  return init;
}

template <ranges::input_range Range>
inline auto accumulate_sum(Range&& range) -> ranges::range_value_t<Range> {
  return accumulate_sum(std::forward<Range>(range), ranges::range_value_t<Range>{});
}

template <class T>
inline auto make_reserved_vector(size_t capacity) {
  auto res = std::vector<T>{};
  res.reserve(capacity);
  return res;
}

namespace details {
template <class Func, class TupleType, size_t... Indices>
inline auto tuple_transform_impl(std::index_sequence<Indices...>, Func&& func, TupleType&& tuple) {
  return std::tuple{std::invoke(func, get<Indices>(std::forward<TupleType>(tuple)))...};
}
} // namespace details

template <class Func, class TupleType>
inline auto tuple_transform(Func&& func, TupleType&& tuple) {
  constexpr auto N = std::tuple_size_v<std::remove_cvref_t<TupleType>>;
  return details::tuple_transform_impl(std::make_index_sequence<N>{}, std::forward<Func>(func),
                                       std::forward<TupleType>(tuple));
}
