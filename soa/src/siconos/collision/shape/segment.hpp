#pragma once

#include "siconos/collision/collision_head.hpp"

namespace siconos::collision::shape {

struct segment : item<> {
  using attributes = gather<
      attribute<"p1", some::vector<some::scalar, some::indice_value<3>>>,
      attribute<"p2", some::vector<some::scalar, some::indice_value<3>>>,
      attribute<"normal", some::vector<some::scalar, some::indice_value<3>>>,
      attribute<"length_square", some::scalar>>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) p1() { return attr<"p1">(*self()); };
    decltype(auto) p2() { return attr<"p2">(*self()); };
    decltype(auto) length_square() { return attr<"length_square">(*self()); };
    decltype(auto) normal() { return attr<"normal">(*self()); };

    void compute_length_square()
    {
      auto v = p2() - p1();
      length_square() =  v * v;
    };

    void compute_normal()
    {
      auto v = p2() - p1();
      normal()[0] = -v[1];
      normal()[1] = v[0];
    }

    decltype(auto) distance(match::vector auto& q)
    {
      auto& x = q(0);
      auto& y = q(1);

      return fabs(attr<"a">(*self()) * x + attr<"b">(*self()) * y +
                  attr<"c">(*self())) *
             attr<"invsqrta2pb2">(*self());
    }

    decltype(auto) projection(match::vector auto& point)
    {
      auto& x = q(0);
      auto& y = q(1);
    }
  };
};
}  // namespace siconos::collision::shape
