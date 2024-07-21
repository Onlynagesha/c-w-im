#pragma once

#include "utils/easylog.h"
#include "utils/utils.h"

namespace details {
template <std::integral T>
struct MinMaxPair {
  using ValueType = T;
  T min;
  T max;
};

template <class T>
struct IsMinMaxPair : std::false_type {};

template <class T>
struct IsMinMaxPair<MinMaxPair<T>> : std::true_type {
  using ValueType = T;
};

template <class Struct, auto Member>
constexpr auto is_radix_sortable_member = //
    member_object_pointer_from<decltype(Member), Struct> &&
    std::is_integral_v<member_type_of_pointer_to_member<decltype(Member)>>;

template <std::default_initializable Struct, auto CurMember, auto... RestMembers,
          class MinMaxPairOrNull = std::nullptr_t>
  requires(is_radix_sortable_member<Struct, CurMember> && (is_radix_sortable_member<Struct, RestMembers> && ...) &&
           (IsMinMaxPair<MinMaxPairOrNull>::value || std::is_same_v<MinMaxPairOrNull, std::nullptr_t>))
inline auto radix_sort_struct_impl(std::span<Struct> to_range, std::span<Struct> from_range,
                                   MinMaxPairOrNull min_max_pair = nullptr) -> void {
  if constexpr (sizeof...(RestMembers) != 0) {
    radix_sort_struct_impl<Struct, RestMembers...>(from_range, to_range); // Swaps
  }
  using ValueType = member_type_of_pointer_to_member<decltype(CurMember)>;
  auto m = to_range.size();
  auto [min, max] = [&]() {
    if constexpr (!std::is_same_v<MinMaxPairOrNull, std::nullptr_t>) {
      return min_max_pair; // O(1) if min, max are known in prior
    } else {
      return ranges::minmax(from_range | views::transform(CurMember)); // O(n) otherwise.
    }
  }();
  auto count = std::vector<size_t>(max - min + 1);

  for (const auto& item : from_range) {
    count[item.*CurMember - min] += 1;
  }
  std::exclusive_scan(count.begin(), count.end(), count.begin(), 0);
  for (const auto& item : from_range) {
    auto& cur_count = count[item.*CurMember - min];
    to_range[cur_count] = item;
    cur_count += 1;
  }
}
} // namespace details

/**
Example: Sort in lexicographical order by the tuple (x, y, z)

struct Point { short x; int y; long z; };
std::vector<Point> points = some_operation();
radix_sort_struct<&Point::x, &Point::y, &Point::z>(points);
 */
template <auto... Members, std::default_initializable Struct>
  requires(std::is_class_v<Struct> && sizeof...(Members) >= 1 &&
           (details::is_radix_sortable_member<Struct, Members> && ...))
inline auto radix_sort_struct(std::span<Struct> range) -> void {
  auto m = range.size();
  auto temp = std::vector<Struct>(m);
  if constexpr (sizeof...(Members) % 2 == 0) {
    // Even: range -> temp -> range
    details::radix_sort_struct_impl<Struct, Members...>(range, temp);
  } else {
    // Odd: range -> temp -> range -> temp, range copy required.
    details::radix_sort_struct_impl<Struct, Members...>(temp, range);
    ranges::copy(temp, range.begin());
  }
}

template <auto... Members, ranges::contiguous_range Range>
  requires(std::is_class_v<ranges::range_value_t<Range>> && sizeof...(Members) >= 1 &&
           !std::is_const_v<std::remove_reference_t<ranges::range_reference_t<Range>>> &&
           (details::is_radix_sortable_member<ranges::range_value_t<Range>, Members> && ...))
inline auto radix_sort_struct(Range&& range) -> void {
  using ValueType = ranges::range_value_t<Range>;
  radix_sort_struct<Members...>(std::span<ValueType>{range});
}

template <auto... Members, std::default_initializable Struct, std::integral ValueType>
  requires(std::is_class_v<Struct> && sizeof...(Members) >= 1 &&
           (details::is_radix_sortable_member<Struct, Members> && ...))
inline auto radix_sort_struct(std::span<Struct> range, ValueType min, ValueType max) -> void {
  auto m = range.size();
  auto temp = std::vector<Struct>(m);
  auto min_max = details::MinMaxPair<ValueType>{min, max};
  if constexpr (sizeof...(Members) % 2 == 0) {
    // Even: range -> temp -> range
    details::radix_sort_struct_impl<Struct, Members...>(range, temp, min_max);
  } else {
    // Odd: range -> temp -> range -> temp, range copy required.
    details::radix_sort_struct_impl<Struct, Members...>(temp, range, min_max);
    ranges::copy(temp, range.begin());
  }
}

template <auto... Members, ranges::contiguous_range Range, std::integral MinMaxValueType>
  requires(std::is_class_v<ranges::range_value_t<Range>> && sizeof...(Members) >= 1 &&
           !std::is_const_v<std::remove_reference_t<ranges::range_reference_t<Range>>> &&
           (details::is_radix_sortable_member<ranges::range_value_t<Range>, Members> && ...))
inline auto radix_sort_struct(Range&& range, MinMaxValueType min, MinMaxValueType max) -> void {
  using RangeValueType = ranges::range_value_t<Range>;
  radix_sort_struct<Members...>(std::span<RangeValueType>{range}, min, max);
}
