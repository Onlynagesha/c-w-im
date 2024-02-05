#pragma once

#include "utils/boost_assert.h"
#include "utils/utils.h"
#include <algorithm>
#include <fmt/color.h>

template <ranges::forward_range Xs>
  requires(std::is_arithmetic_v<ranges::range_value_t<Xs>>)
inline auto make_histogram(Xs&& values, size_t display_width, size_t display_height) -> std::string {
  BOOST_ASSERT_MSG(!ranges::empty(values), "# of values must be non-zero.");
  BOOST_ASSERT_MSG(display_width > 0, "# of rows must be a positive integer.");
  BOOST_ASSERT_MSG(display_height > 0, "# of columns must be a positive integer.");

  // U+2581 ~ U+2588: each incrementing 1/8 height
  // Note: all the string literals below are in UTF-8 encoding.
  constexpr auto BLOCKS = std::array{"▁", "▂", "▃", "▄", "▅", "▆", "▇", "█"};

  auto [x_min, x_max] = ranges::minmax(values);
  auto x_seg_length = 1.0 * (x_max - x_min) / display_width;

  auto heights = std::vector<double>(display_width, 0.0);
  for (auto val : values) {
    auto idx = std::clamp(static_cast<size_t>((val - x_min) / x_seg_length), 0zu, heights.size() - 1);
    heights[idx] += 1.0;
  }
  auto [height_min, height_max] = ranges::minmax(heights);
  auto to_display_height = [&](auto h) {
    auto h_disp = size(BLOCKS) * display_height * h / height_max;
    return static_cast<size_t>(std::ceil(h_disp)); // Rounds up, and converts to UNSIGNED integer
  };
  auto display_heights = [&] {
    auto view =  heights | views::transform(to_display_height);
    return std::vector(view.begin(), view.end());
  }();

  auto horizontal_line = [&] {
    auto view = views::repeat("─"sv, display_width) | views::join;
    return std::string(view.begin(), view.end());
  }();
  auto res = "┌" + horizontal_line + "┐";
  for (auto base_line : range(display_height) | views::transform(LAMBDA_1(size(BLOCKS) * _1)) | views::reverse) {
    res += "\n│";
    for (auto height : display_heights) {
      if (height > base_line) {
        res += BLOCKS[std::min(height - base_line, size(BLOCKS)) - 1];
      } else {
        res += " ";
      }
    }
    res += "│";
  }

  auto n = 0.0;
  auto sum = 0.0;
  for (auto val : values) {
    n += 1.0; // In case of Xs is not sized range
    sum += val;
  }
  auto mean = sum / n;
  auto sqr_mean_diff = views::transform([mean](auto x) { return (x - mean) * (x - mean); });
  auto std = std::sqrt(accumulate_sum(values | sqr_mean_diff) / n);

  constexpr auto final_pattern = "\n└{}┘\n[x in range [{}, {}], y in range [{}, {}]], "
                                 "# = {}, mean = {:.3f}, std = {:.3f}";
  res += fmt::format(final_pattern, horizontal_line, x_min, x_max, height_min, height_max, n, mean, std);
  return res;
}
