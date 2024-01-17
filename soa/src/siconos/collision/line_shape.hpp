#pragma once

#include "siconos/collision/collision.hpp"

namespace siconos::collision {

struct line_shape : item<> {
  // a*x + b*y + c = 0
  using attributes =
      gather<attribute<"a", some::scalar>, attribute<"b", some::scalar>,
             attribute<"c", some::scalar>,
             attribute<"invsqrta2pb2", some::scalar>>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) distance(match::vector auto& q)
    {
      auto& x = q(0);
      auto& y = q(1);

      return fabs(attr<"a">(*self()) * x + attr<"b">(*self()) * y +
                  attr<"c">(*self())) *
             attr<"invsqrta2pb2">(*self());
    }
  };
};
}  // namespace siconos::collision
