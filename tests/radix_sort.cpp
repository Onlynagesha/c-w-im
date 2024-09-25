#define BOOST_TEST_MODULE "Radix sort"

#include "radix_sort.h"
#include <boost/test/unit_test.hpp>
#include <fmt/ranges.h>
#include <nwgraph/util/timer.hpp>

struct Point1D {
  uint8_t x;

  constexpr auto operator<=>(const Point1D& rhs) const -> std::strong_ordering = default;
  constexpr auto operator==(const Point1D& rhs) const -> bool = default;

  friend auto operator<<(std::ostream& out, const Point1D& p) -> std::ostream& {
    out << fmt::format("({})", p.x);
    return out;
  }
};

BOOST_AUTO_TEST_CASE(point1D) {
  auto points = std::vector<Point1D>{
      {1}, {3}, {6}, {9}, {static_cast<uint8_t>(-1)}, {1}, {2}, {5}, {7}, {static_cast<uint8_t>(-2)},
  };
  auto points_copy = points;

  auto timer = nw::util::us_timer();
  timer.start();
  radix_sort_struct<&Point1D::x>(points);
  timer.stop();
  std::cout << "Point1D list after radix sorting:" << std::endl;
  for (auto p : points) {
    std::cout << "\t" << p << std::endl;
  }

  ranges::sort(points_copy);
  BOOST_CHECK_EQUAL_COLLECTIONS(points.begin(), points.end(), points_copy.begin(), points_copy.end());
  std::cout << fmt::format("Point1D test case done. Time used: {:.3f} us.", timer.elapsed()) << std::endl;
}

struct Point2D {
  short x;
  int y;

  constexpr auto operator<=>(const Point2D& rhs) const -> std::strong_ordering = default;
  constexpr auto operator==(const Point2D& rhs) const -> bool = default;

  friend auto operator<<(std::ostream& out, const Point2D& p) -> std::ostream& {
    out << fmt::format("({}, {})", p.x, p.y);
    return out;
  }
};

BOOST_AUTO_TEST_CASE(point2D) {
  auto points = std::vector<Point2D>{
      {0, 1}, {2, -5}, {3, 6}, {-3, 1}, {6, 2}, {-1, 2}, {7, 3}, {0, -8}, {1, 3}, {5, -3},
      {6, 2}, {0, -1}, {1, 4}, {3, -5}, {3, 0}, {2, -8}, {1, 1}, {-3, 5}, {2, 9}, {-2, 5},
  };
  auto points_copy = points;

  auto timer = nw::util::us_timer();
  timer.start();
  radix_sort_struct<&Point2D::x, &Point2D::y>(points);
  timer.stop();
  std::cout << "Point2D list after radix sorting:" << std::endl;
  for (auto p : points) {
    std::cout << "\t" << p << std::endl;
  }

  ranges::sort(points_copy);
  BOOST_CHECK_EQUAL_COLLECTIONS(points.begin(), points.end(), points_copy.begin(), points_copy.end());
  std::cout << fmt::format("Point2D test case done. Time used: {:.3f} us.", timer.elapsed()) << std::endl;
}

struct Point3D {
  int16_t x;
  int32_t y;
  uint64_t z;

  constexpr auto operator<=>(const Point3D& rhs) const -> std::strong_ordering = default;
  constexpr auto operator==(const Point3D& rhs) const -> bool = default;

  friend auto operator<<(std::ostream& out, const Point3D& p) -> std::ostream& {
    out << fmt::format("({}, {}, {})", p.x, p.y, p.z);
    return out;
  }
};

BOOST_AUTO_TEST_CASE(point3D) {
  auto points = std::vector<Point3D>{
      {0, -1, 3}, {3, 2, 1}, {-3, 4, 2}, {1, 3, 5}, {2, -3, 5}, {3, 2, 1}, {-2, 5, 9}, {2, 1, 4}, {3, 4, 1},
      {4, -2, 9}, {3, 1, 9}, {8, -3, 1}, {2, 5, 2}, {4, -2, 3}, {5, 3, 3}, {-1, 0, 3}, {0, 1, 8}, {2, 3, 1},
      {-1, 9, 3}, {1, 2, 0}, {3, -1, 4}, {5, 2, 1}, {-3, 5, 1}, {2, 5, 4}, {5, -1, 0},
  };
  auto points_copy = points;

  auto timer = nw::util::us_timer();
  timer.start();
  radix_sort_struct<&Point3D::x, &Point3D::y, &Point3D::z>(points);
  timer.stop();
  std::cout << "Point3D list after radix sorting:" << std::endl;
  for (auto p : points) {
    std::cout << "\t" << p << std::endl;
  }

  ranges::sort(points_copy);
  BOOST_CHECK_EQUAL_COLLECTIONS(points.begin(), points.end(), points_copy.begin(), points_copy.end());
  std::cout << fmt::format("Point3D test case done. Time used: {:.3f} us.", timer.elapsed()) << std::endl;
}
