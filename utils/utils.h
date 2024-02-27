#pragma once

#include "utils/boost_assert.h"
#include <concepts>
#include <fmt/format.h>
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

#define FILTER_VIEW(...) std::views::filter(LAMBDA_1(__VA_ARGS__))
#define TRANSFORM_VIEW(...) std::views::transform(LAMBDA_1(__VA_ARGS__))

#define DUMP_ARRAY(arr) dump_array(#arr, arr)
#define DUMP_INDEX_ARRAY(arr) dump_index_array(#arr, arr)

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

template <class T, size_t N>
using Array1D = std::array<T, N>;

template <class T, size_t Rows, size_t Columns>
using Array2D = std::array<std::array<T, Columns>, Rows>;

// Alternative of std::forward_like: https://en.cppreference.com/w/cpp/utility/forward_like
template <class T, class U>
inline constexpr auto&& forward_like(U&& x) noexcept {
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

namespace details {
template <auto K, class T>
struct IntegerSequenceOffsetHelper {};

template <auto K, class T, T... Values>
struct IntegerSequenceOffsetHelper<K, std::integer_sequence<T, Values...>> {
  using type = std::integer_sequence<T, (static_cast<T>(K) + Values)...>;
};

template <class T>
struct MemberObjectPointerHelper : std::false_type {};

template <class ClassType, class MemberType>
struct MemberObjectPointerHelper<MemberType ClassType::*> : std::true_type {
  using class_type = ClassType;
  using member_type = MemberType;
};
} // namespace details

template <size_t N, size_t Offset>
using make_index_sequence_by_offset = details::IntegerSequenceOffsetHelper<Offset, std::make_index_sequence<N>>::type;

template <class T, class ClassType>
concept member_object_pointer_from =
    details::MemberObjectPointerHelper<T>::value &&
    std::is_same_v<typename details::MemberObjectPointerHelper<T>::class_type, ClassType>;

template <class T, class ClassType, class MemberType>
concept member_object_pointer_from_to =
    details::MemberObjectPointerHelper<T>::value &&
    std::is_same_v<typename details::MemberObjectPointerHelper<T>::class_type, ClassType> &&
    std::is_same_v<typename details::MemberObjectPointerHelper<T>::member_type, MemberType>;

#define UTILS_RANGE_CATEGORIES(F) \
  F(input)                        \
  F(forward)                      \
  F(bidirectional)                \
  F(random_access)                \
  F(contiguous)

#define UTILS_DEFINE_CONCEPT_RANGE_OF(category)                                                                       \
  template <class T, class ValueType>                                                                                 \
  concept category##_range_of =                                                                                       \
      ranges::category##_range<T> && std::is_convertible_v<ranges::range_value_t<std::remove_cvref_t<T>>, ValueType>; \
  template <class T, class ValueType>                                                                                 \
  concept category##_range_of_exactly =                                                                               \
      ranges::category##_range<T> && std::is_same_v<ranges::range_value_t<std::remove_cvref_t<T>>, ValueType>;

UTILS_RANGE_CATEGORIES(UTILS_DEFINE_CONCEPT_RANGE_OF)

// WARNING: No thread safety.
#ifdef DEBUG_FIXED_RANDOM_SEED
inline auto rand_engine = std::minstd_rand{1u};
#else
inline auto rand_engine = std::minstd_rand{std::random_device{}()};
#endif

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
  BOOST_ASSERT_MSG(0.0 <= p && p <= 1.0, "p must be in the range [0, 1].");
  return rand_float() < p;
}

inline auto rand_index(size_t n) -> size_t {
  BOOST_ASSERT_MSG(n > 0, "n must be a positive integer.");
  auto dist = std::uniform_int_distribution<size_t>{0, n - 1};
  return dist(rand_engine);
}

template <class Range>
  requires(ranges::random_access_range<std::remove_cvref_t<Range>> && ranges::sized_range<std::remove_cvref_t<Range>>)
inline auto rand_element(Range&& range) -> ranges::range_value_t<std::remove_cvref_t<Range>> {
  return range[rand_index(ranges::size(range))];
}

template <class T>
inline constexpr auto remove_const(const T& value) -> T& {
  return const_cast<T&>(value);
}

template <std::integral IntType>
inline constexpr auto range(IntType n) {
  return views::iota(static_cast<IntType>(0), n);
}

template <std::integral IntType>
inline constexpr auto range(IntType first, IntType last) {
  return views::iota(first, last);
}

template <std::input_or_output_iterator Iter, std::sentinel_for<Iter> Sentinel>
inline constexpr auto range(Iter first, Sentinel last) {
  return ranges::subrange{first, last};
}

template <ranges::input_range Range>
  requires(std::unsigned_integral<ranges::range_value_t<Range>>)
inline constexpr auto make_bitset_from_indices(Range&& indices) -> ranges::range_value_t<Range> {
  using ResultType = ranges::range_value_t<Range>;
  auto res = ResultType{0};
  for (auto i : indices) {
    res |= (ResultType{1} << i);
  }
  return res;
}

template <std::unsigned_integral T>
inline constexpr auto make_bitset_from_indices(std::initializer_list<T> indices) -> T {
  return make_bitset_from_indices(std::span{indices});
}

template <std::unsigned_integral T>
inline constexpr auto bitset_contains_index(T index, T bitset) -> bool {
  return (bitset & (T{1} << index)) != 0;
}

template <std::unsigned_integral T>
inline constexpr auto indices_in_bitset(T N, T bitset) {
  return range(N) | std::views::filter([bitset](T index) { return bitset_contains_index(index, bitset); });
}

template <std::unsigned_integral T>
inline constexpr auto indices_in_bitset(T bitset) {
  return indices_in_bitset(std::bit_width(bitset), bitset);
}

template <class T, ranges::input_range Range>
  requires(std::is_convertible_v<ranges::range_value_t<Range>, T>)
inline constexpr auto accumulate_sum(Range&& range, T init) -> T {
  for (auto&& elem : range) {
    init = init + forward_like<Range>(elem);
  }
  return init;
}

template <ranges::input_range Range>
inline constexpr auto accumulate_sum(Range&& range) -> ranges::range_value_t<Range> {
  return accumulate_sum(range, ranges::range_value_t<Range>{});
}

template <class T, ranges::input_range Range>
  requires(std::is_convertible_v<ranges::range_value_t<Range>, T>)
inline constexpr auto accumulate_product(Range&& range, T init) -> T {
  for (auto&& elem : range) {
    init = init * forward_like<Range>(elem);
  }
  return init;
}

template <ranges::input_range Range>
inline constexpr auto accumulate_product(Range&& range) -> ranges::range_value_t<Range> {
  return accumulate_product(range, ranges::range_value_t<Range>{1});
}

template <std::floating_point... Args>
inline constexpr auto at_least_1_probability(Args... p_args) -> std::common_type_t<Args...> {
  BOOST_ASSERT_MSG(((0.0 <= p_args && p_args <= 1.0) && ...), "p must be in the range [0, 1].");
  using ResultType = std::common_type_t<Args...>;
  auto prod = (static_cast<ResultType>(1.0) * ... * (static_cast<Args>(1.0) - p_args));
  return static_cast<ResultType>(1.0) - prod;
}

template <ranges::input_range Range>
  requires(std::is_floating_point_v<ranges::range_value_t<Range>>)
inline constexpr auto at_least_1_probability_of_range(Range&& range) -> ranges::range_value_t<Range> {
  using ResultType = ranges::range_value_t<Range>;
  auto prod = static_cast<ResultType>(1.0);
  for (auto p : range) {
    BOOST_ASSERT_MSG(0.0 <= p && p <= 1.0, "p must be in the range [0, 1].");
    prod *= (static_cast<ResultType>(1.0) - p);
  }
  return static_cast<ResultType>(1.0) - prod;
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

template <class T>
auto dump_span(std::string_view name, std::span<const T> values) -> std::string {
  auto res = fmt::format("Elements of array '{}' (size = {}):", name, values.size());
  for (const auto& [i, x] : views::enumerate(values)) {
    res += fmt::format("\n\t{}[{}] = {}", name, i, x);
  }
  return res;
}

template <std::unsigned_integral T>
auto dump_index_span(std::string_view name, std::span<const T> values) -> std::string {
  auto res = fmt::format("Elements of index array '{}' (size = {}):", name, values.size());
  for (const auto& [i, x] : views::enumerate(values)) {
    if (x != static_cast<T>(-1)) {
      res += fmt::format("\n\t{}[{}] = {}", name, i, x);
    } else {
      res += fmt::format("\n\t{}[{}] = NULL", name, i);
    }
  }
  return res;
}

template <ranges::random_access_range Range>
auto dump_array(std::string_view name, const Range& arr) -> std::string {
  using ValueType = ranges::range_value_t<Range>;
  return dump_span(name, std::span<const ValueType>{arr.begin(), arr.end()});
}

template <ranges::random_access_range Range>
auto dump_index_array(std::string_view name, const Range& arr) -> std::string {
  using ValueType = ranges::range_value_t<Range>;
  return dump_index_span(name, std::span<const ValueType>{arr.begin(), arr.end()});
}
