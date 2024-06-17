#pragma once

#include "siconos/collision/collision_head.hpp"

namespace siconos::collision::shape {
struct disk : item<> {
  using attributes = gather<attribute<"radius", some::scalar>>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) radius() { return attr<"radius">(*self()); };
  };
};
}  // namespace siconos::collision
