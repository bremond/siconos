#pragma once

#include "siconos/utils/pattern.hpp"

namespace siconos {

struct lagrangian_r : item<> {

  
  using attributes = types::attributes<>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    template <match::any_full_handle DynamicalSystem,
              match::any_full_handle Interaction>
    void compute_input(double time, DynamicalSystem ds, Interaction inter,
                       auto level)
    {
      auto &lambda = inter.lambda()[level][0];
      prod(trans(inter.h_matrix()[0]), lambda,
           ds.property(symbol<"p0">{})[0]);
      // link to ds variables => check graph
    };

    template <match::any_full_handle DynamicalSystem,
              match::any_full_handle Interaction>
    void compute_output(double time, DynamicalSystem ds, Interaction inter,
                        auto level)
    {
      auto &y = inter.y()[1];
      auto &H = inter.h_matrix()[0];

      auto &x = [&level, &ds]() -> decltype(auto) {
        if constexpr (level == 1) {
          return ds.velocity();
        }
        else {
          return ds.q();
          ;
        }
      }();

      prod(H, x, y);
    };
  };
};
}  // namespace siconos
