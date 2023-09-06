#pragma once

#include "siconos/storage/storage.hpp"
#include "siconos/utils/pattern.hpp"

namespace siconos::model {

struct lagrangian_r : item<> {
  using attributes = types::attributes<>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;
  };
};

struct lagrangian_tir
    : item<description<"A lagrangian time invariant relation">> {
  using dof = some::indice_parameter<"dof">;

//  struct h_matrix : some::matrix<some::scalar, nslaw_size, dof>,
  //                   access<h_matrix> {};

  struct b : some::scalar, access<b> {};

  using attributes = gather<b>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    // decltype(auto) h_matrix() { return Handle::type::h_matrix::at(*self());}

    decltype(auto) b() { return Handle::type::b::at(*self());}
  };
};

}  // namespace siconos::model
