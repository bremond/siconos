#pragma once

#include "siconos/collision/collision.hpp"

namespace siconos::collision {
struct disk_shape : item<> {
  using attributes = gather<attribute<"radius", some::scalar>>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) radius() { return attr<"radius">(*self()); };
  };
};
}  // namespace siconos::collision
