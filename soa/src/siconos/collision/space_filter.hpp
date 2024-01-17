#pragma once

#include "siconos/collision/collision.hpp"
#include "siconos/collision/disk_shape.hpp"
#include "siconos/collision/line_shape.hpp"

namespace siconos::collision {
struct space_filter : item<> {
  using attributes =
      gather<attribute<"lines", some::unbounded_collection<line_shape>>,
             attribute<"disks", some::unbounded_collection<disk_shape>>>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    decltype(auto) insert_line(auto a, auto b, auto c)
    {
      attr<"lines">(self()).push_back({a, b, c});
    }
  };
}
}  // namespace siconos::collision
