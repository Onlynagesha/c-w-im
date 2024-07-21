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

#define CHUNK_BY_VIEW(...) std::views::chunk_by(LAMBDA_2(__VA_ARGS__))
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

// Gets std::index_sequence<Offset, Offset+1 ... Offset+(N-1)>
template <size_t N, size_t Offset>
using make_index_sequence_by_offset = details::IntegerSequenceOffsetHelper<Offset, std::make_index_sequence<N>>::type;

// Checks whether T is a pointer to member from given ClassType
// e.g. struct Point { int x; float y; double z; };
//      static_assert( member_object_pointer_from<&Point::x, Point> );
template <class T, class ClassType>
concept member_object_pointer_from =
    details::MemberObjectPointerHelper<T>::value &&
    std::is_same_v<typename details::MemberObjectPointerHelper<T>::class_type, ClassType>;

// Checks whether T is a pointer to member of EXACTLY type MemberType from given ClassType
// e.g. struct Point { int x; float y; double z; };
//      static_assert( member_object_pointer_from<&Point::y, Point, float> );
template <class T, class ClassType, class MemberType>
concept member_object_pointer_from_to =
    details::MemberObjectPointerHelper<T>::value &&
    std::is_same_v<typename details::MemberObjectPointerHelper<T>::class_type, ClassType> &&
    std::is_same_v<typename details::MemberObjectPointerHelper<T>::member_type, MemberType>;

template <class T>
using class_type_of_pointer_to_member = typename details::MemberObjectPointerHelper<T>::class_type;

template <class T>
using member_type_of_pointer_to_member = typename details::MemberObjectPointerHelper<T>::member_type;

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

// Examples:
// static_assert( forward_range_of<std::forward_list<float>, double> );     // Implicit conversion is allowed
// static_assert( bidirectional_range_of_exactly<std::set<float>, float> ); // Exact type match is required
// static_assert( random_acess_range_of<std::deque<unsigned>, size_t> );
// static_assert( contiguous_range_of_exactly<std::array<unsigned, 32>, unsigned> );
UTILS_RANGE_CATEGORIES(UTILS_DEFINE_CONCEPT_RANGE_OF)

// WARNING: No thread safety.
// Linear algorithm is used for the sake of performance.
#ifdef DEBUG_FIXED_RANDOM_SEED
inline auto rand_engine = std::minstd_rand{1u};
#else
inline auto rand_engine = std::minstd_rand{std::random_device{}()};
#endif

// Generates a random float in range [0, 1).
inline auto rand_float() -> float {
  constexpr auto ieee754_float_fraction_bits = 23u;
  constexpr auto ieee754_float_fraction_mask = (1u << ieee754_float_fraction_bits) - 1;
  constexpr auto ieee754_float_non_fraction_part = 127u << ieee754_float_fraction_bits;

  auto u = rand_engine();
  // res as float32 = 1.??? x 2^0 = random value in range [1, 2)
  auto res = (u & ieee754_float_fraction_mask) | ieee754_float_non_fraction_part;
  return reinterpret_cast<float&>(res) - 1;
}

// Generates true with probability p, false with 1-p.
inline auto rand_bool(float p) -> bool {
  BOOST_ASSERT_MSG(0.0 <= p && p <= 1.0, "p must be in the range [0, 1].");
  return rand_float() < p;
}

// Generates a random index uniformly in range [0, n).
inline auto rand_index(size_t n) -> size_t {
  BOOST_ASSERT_MSG(n > 0, "n must be a positive integer.");
  auto dist = std::uniform_int_distribution<size_t>{0, n - 1};
  return dist(rand_engine);
}

// Takes a random element by uniform distribution from the given range.
// Returns a COPY of the randomly selected element.
template <class Range>
  requires(ranges::random_access_range<Range> && ranges::sized_range<Range>)
inline auto rand_element(Range&& range) -> ranges::range_value_t<Range> {
  return range[rand_index(ranges::size(range))];
}

struct as_non_const_t {
  template <class T>
  constexpr auto operator()(const T& value) const -> T& {
    return const_cast<T&>(value);
  }
};
// Equivalent to const_cast<T&>(x) where x is a const reference
constexpr auto as_non_const = as_non_const_t{};

struct to_signed_t {
  template <std::integral T>
  constexpr auto operator()(T value) const {
    using S = std::make_signed_t<T>;
    return static_cast<S>(value);
  }
};
// Equivalent to static_cast<S>(u) where u is an unsigned integer and S is its corresponding signed type.
constexpr auto to_signed = to_signed_t{};

struct to_unsigned_t {
  template <std::integral T>
  constexpr auto operator()(T value) const {
    using U = std::make_unsigned_t<T>;
    return static_cast<U>(value);
  }
};
// Equivalent to static_cast<U>(s) where s is a signed integer and U is its corresponding unsigned type.
constexpr auto to_unsigned = to_unsigned_t{};

// Equivalent to range(n) in Python, or std::views::iota(0, n) in C++.
template <std::integral IntType>
inline constexpr auto range(IntType n) {
  return views::iota(static_cast<IntType>(0), n);
}

// Equivalent to range(first, last) in Python, or std::views::iota(first, last) in C++.
// first and last must be of exactly the same integer type.
template <std::integral IntType>
inline constexpr auto range(IntType first, IntType last) {
  return views::iota(first, last);
}

// Equivalent to std::ranges::subrange(first, last) in C++.
template <std::input_or_output_iterator Iter, std::sentinel_for<Iter> Sentinel>
inline constexpr auto range(Iter first, Sentinel last) {
  return ranges::subrange{first, last};
}

// Transforms index sequence to bitset.
// e.g. indices = [0, 2, 3, 5, 8], then result = 0b1'0010'1101
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

// Transforms index sequence (given as std::initializer_list) to bitset.
// e.g. indices = [0, 2, 3, 5, 8], then result = 0b1'0010'1101
template <std::unsigned_integral T>
inline constexpr auto make_bitset_from_indices(std::initializer_list<T> indices) -> T {
  return make_bitset_from_indices(std::span{indices});
}

// Tests whether bitset contains the given index.
template <std::unsigned_integral T>
inline constexpr auto bitset_contains_index(T index, T bitset) -> bool {
  return (bitset & (T{1} << index)) != 0;
}

// Generates a range of indices contained in bitset.
// e.g. N = 12, bitset = 0b0001'0010'1101, then the result range is [0, 2, 3, 5, 8].
template <std::unsigned_integral T>
inline constexpr auto indices_in_bitset(T N, T bitset) {
  return range(N) | std::views::filter([bitset](T index) { return bitset_contains_index(index, bitset); });
}

// Generates a range of indices contained in bitset, with N deduced automatically.
template <std::unsigned_integral T>
inline constexpr auto indices_in_bitset(T bitset) {
  return indices_in_bitset(std::bit_width(bitset), bitset);
}

// Equivalent to std::accumulate(begin(range), end(range), init)
template <class T, ranges::input_range Range>
  requires(std::is_convertible_v<ranges::range_value_t<Range>, T>)
inline constexpr auto accumulate_sum(Range&& range, T init) -> T {
  for (auto&& elem : range) {
    init = init + forward_like<Range>(elem);
  }
  return init;
}

// Equivalent to std::accumulate(begin(range), end(range), T{}) where T is the value type of range.
// For arithmetic types, init = 0.
template <ranges::input_range Range>
inline constexpr auto accumulate_sum(Range&& range) -> ranges::range_value_t<Range> {
  return accumulate_sum(range, ranges::range_value_t<Range>{});
}

// Equivelent to std::accumulate(begin(range), end(range), init, std::multiplies{}).
template <class T, ranges::input_range Range>
  requires(std::is_convertible_v<ranges::range_value_t<Range>, T>)
inline constexpr auto accumulate_product(Range&& range, T init) -> T {
  for (auto&& elem : range) {
    init = init * forward_like<Range>(elem);
  }
  return init;
}

// Equivalent to std::accumulate(begin(range), end(range), 1, std::multiplies{}).
template <ranges::input_range Range>
  requires(std::is_arithmetic_v<ranges::range_value_t<Range>>)
inline constexpr auto accumulate_product(Range&& range) -> ranges::range_value_t<Range> {
  return accumulate_product(range, ranges::range_value_t<Range>{1});
}

// Given N independent events of probability p1 ... pN,
// Gets the probability that at least 1 of the events happens,
// where N = sizeof...(Args), p1 ... pN is expanded from p_args.
template <std::floating_point... Args>
inline constexpr auto at_least_1_probability(Args... p_args) -> std::common_type_t<Args...> {
  BOOST_ASSERT_MSG(((0.0 <= p_args && p_args <= 1.0) && ...), "p must be in the range [0, 1].");
  using ResultType = std::common_type_t<Args...>;
  auto prod = (static_cast<ResultType>(1.0) * ... * (static_cast<Args>(1.0) - p_args));
  return static_cast<ResultType>(1.0) - prod;
}

// Given N independent events of probability p1 ... pN,
// Gets the probability that at least 1 of the events happens,
// where N = size(range), range = [p1 ... pN].
template <ranges::input_range Range>
  requires(std::is_floating_point_v<ranges::range_value_t<Range>>)
inline constexpr auto at_least_1_probability_r(Range&& range) -> ranges::range_value_t<Range> {
  using ResultType = ranges::range_value_t<Range>;
  auto prod = static_cast<ResultType>(1.0);
  for (auto p : range) {
    BOOST_ASSERT_MSG(0.0 <= p && p <= 1.0, "p must be in the range [0, 1].");
    prod *= (static_cast<ResultType>(1.0) - p);
  }
  return static_cast<ResultType>(1.0) - prod;
}

//   auto vec = make_reserved_vector<int>(256);
// Equivalent to the following:
//   auto vec = std::vector<int>();
//   vec.reserve(256);
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

// Applies transform to a tuple.
// Example:
// auto func = [](auto x) { return x + x; };
// auto [i, f, s] = tuple_transform(func, std::tuple{21, 1.5, "abc"s});
// assert(i == 42 && f == 3.0 && s == "abcabc");
template <class Func, class TupleType>
inline auto tuple_transform(Func&& func, TupleType&& tuple) {
  constexpr auto N = std::tuple_size_v<std::remove_cvref_t<TupleType>>;
  return details::tuple_transform_impl(std::make_index_sequence<N>{}, std::forward<Func>(func),
                                       std::forward<TupleType>(tuple));
}

// Dumps value range
template <class T>
inline auto dump_span(std::string_view name, std::span<const T> values) -> std::string {
  auto res = fmt::format("Elements of array '{}' (size = {}):", name, values.size());
  for (const auto& [i, x] : views::enumerate(values)) {
    res += fmt::format("\n\t{}[{}] = {}", name, i, x);
  }
  return res;
}

// Dumps index range (with -1 replaced by "NULL")
template <std::unsigned_integral T>
inline auto dump_index_span(std::string_view name, std::span<const T> values) -> std::string {
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

// Dumps array (see above)
template <ranges::contiguous_range Range>
inline auto dump_array(std::string_view name, const Range& arr) -> std::string {
  using ValueType = ranges::range_value_t<Range>;
  return dump_span(name, std::span<const ValueType>{arr.begin(), arr.end()});
}

// Dumps index array (see above, -1 replaced by "NULL")
template <ranges::contiguous_range Range>
inline auto dump_index_array(std::string_view name, const Range& arr) -> std::string {
  using ValueType = ranges::range_value_t<Range>;
  return dump_index_span(name, std::span<const ValueType>{arr.begin(), arr.end()});
}

// Converts size bytes to its string representation.
// e.g.     123'456'789 => "117.738 Mebibytes"
//      123'456'789'012 => "114.978 Gibibytes"
inline auto size_bytes_to_memory_str(size_t size_bytes) -> std::string {
  constexpr auto UNITS = std::array{"Bytes", "KibiBytes", "Mebibytes", "Gibibytes"};
  auto value = static_cast<double>(size_bytes);
  auto unit_index = 0;
  while (value >= 1024.0 && unit_index + 1 < UNITS.size()) {
    value /= 1024.0;
    unit_index += 1;
  }
  return fmt::format("{:.3f} {}", value, UNITS[unit_index]);
}

// Estimates memory allocation of std::vector<T>.
// Nested std::vector is supported. e.g. std::vector<std::vector<int>>
// Yet nested range of another type is not. e.g. std::vector<std::deque<int>>
template <class T>
inline auto estimated_memory_usage(const std::vector<T>& vector) -> size_t {
  auto res = sizeof(T) * vector.capacity();
  if constexpr (is_std_vector_v<T>) {
    for (const auto& inner : vector) {
      res += estimated_memory_usage(inner);
    }
  } else {
    static_assert(!ranges::range<T>, "Estimation of nested range only works for std::vector.");
  }
  return res;
}
