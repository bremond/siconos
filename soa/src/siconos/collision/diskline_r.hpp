#pragma once

#include "siconos/collision/collision.hpp"

namespace siconos::collision {

  struct diskline_r : item<>, model::relation1, model::any_lagrangian_relation {
  using attributes = gather<attribute<"disk", some::item_ref<disk_shape>>,
                            attribute<"line", some::item_ref<line_shape>>>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) r()
    {
      return handle(self()->data(), attr<"disk">(*self())).radius();
    };
    decltype(auto) line()
    {
      return handle(self()->data(), attr<"line">(*self()));
    };

    decltype(auto) compute_h(match::vector auto& q)
    {
      return line().distance(q);
    }

    decltype(auto) compute_jachq(auto step, auto& q, auto& h_matrix1)
    {
      auto& x = q(0);
      auto& y = q(1);
      auto& a = attr<"a">(line());
      auto& b = attr<"b">(line());
      auto& c = attr<"c">(line());
      auto& invsqrta2pb2 = attr<"invsqrta2pb2">(line());

      auto d1 = a * x + b * y + c;
      auto sd1 = copysign(1, d1);

      h_matrix1(0, 0) = a * sd1 * invsqrta2pb2;
      h_matrix1(1, 0) = -b * sd1 * invsqrta2pb2;
      h_matrix1(0, 1) = b * sd1 * invsqrta2pb2;
      h_matrix1(1, 1) = a * sd1 * invsqrta2pb2;
      h_matrix1(0, 2) = 0;
      h_matrix1(1, 2) = -r();
    }
  };
};
}
