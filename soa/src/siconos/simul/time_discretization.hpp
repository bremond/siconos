#pragma once

#include "siconos/siconos.hpp"
#include "siconos/utils/pattern.hpp"

namespace siconos::simul {

template <typename... Params>
struct time_discretization : item<> {
  struct t0 : some::scalar, access<t0> {};
  struct tmax : some::scalar, access<tmax> {};
  struct h : some::scalar, access<h> {};
  struct step : some::indice, access<step> {};
  using attributes = gather<h, t0, tmax, step>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) t0() { return Handle ::type ::t0 ::at(*self()); };
    decltype(auto) tmax() { return Handle ::type ::tmax ::at(*self()); };
    decltype(auto) h() { return Handle ::type ::h ::at(*self()); };
    decltype(auto) step() { return Handle ::type ::step ::at(*self()); };
  };
};
}  // namespace siconos
