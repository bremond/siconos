#pragma once

#include "siconos/collision/collision.hpp"

namespace siconos::collision {

  struct diskdisk_r : item<>, model::relation2, model::any_lagrangian_relation {
  using dof = some::indice_parameter<"dof">;

  using attributes = gather<attribute<"disk1", some::item_ref<disk_shape>>,
                            attribute<"disk2", some::item_ref<disk_shape>>>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) r1()
    {
      return handle(self()->data(), attr<"disk1">(*self())).radius();
    };
    decltype(auto) r2()
    {
      return handle(self()->data(), attr<"disk2">(*self())).radius();
    };

    template <match::vector V, match::scalar S>
    S compute_h(V& q1, V& q2)
    {
      auto dx = q2[0] - q1[0];
      auto dy = q2[1] - q2[1];

      return algebra::hypot(dx, dy) - (r1() + r2());
    }

    template <typename S, typename V, typename M>
    decltype(auto) compute_jachq(S step, V& q1, V& q2, M& h_matrix1,
                                 M& h_matrix2)
    {
      auto x1 = q1(0);
      auto y1 = q1(1);
      auto x2 = q2(0);
      auto y2 = q2(1);

      auto dx = x2 - x1;
      auto dy = y2 - y1;

      auto d = algebra::hypot(dx, dy);

      auto dxsd = dx / d;
      auto dysd = dy / d;

      auto& g1 = h_matrix1;
      auto& g2 = h_matrix2;

      g1(0, 0) = -dxsd;
      g1(1, 0) = dysd;
      g1(0, 1) = -dysd;
      g1(1, 1) = -dxsd;
      g1(0, 2) = 0.;
      g1(1, 2) = -r1();

      g2(0, 0) = dxsd;
      g2(1, 0) = -dysd;
      g2(0, 1) = dysd;
      g2(1, 1) = dxsd;
      g2(0, 2) = 0.;
      g2(1, 2) = -r2();
    }

    auto methods()
    {
      using env_t = decltype(self()->env());
      using scalar = typename env_t::scalar;
      using params_t = typename env_t::params;

      constexpr auto dof =
          ground::get<pattern::param<"dof">>(params_t{}).value;

      using vector = typename env_t::template vector<scalar, dof>;

      return collect(
          method("compute_h", &interface<Handle>::compute_h<vector, scalar>),
          method("r1", &interface<Handle>::r1),
          method("r2", &interface<Handle>::r2)
          //          ground::make_key_value(
          //              "compute_jachq"_s,
          //              &interface<Handle>::compute_jachq<scalar, vector,
          //              matrix>
      );
    }
  };

  //             method<"compute_jachq", &interface<H>::compute_jachq>>;
  // cf
  // https://stackoverflow.com/questions/47906426/type-of-member-functions-arguments
};
}  // namespace siconos::collision
