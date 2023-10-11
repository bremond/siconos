#pragma once

#include "siconos/model/model.hpp"

namespace siconos::model {

struct lagrangian_ds
    : item<description<"A lagrangian dynamical system [...]">> {
  using dof = some::indice_parameter<"dof">;

  struct q : some::vector<some::scalar, dof>, access<q> {};

  struct velocity : some::vector<some::scalar, dof>, access<velocity> {};

  // should not be an attribute
  struct fext : some::vector<some::scalar, dof>,  // some::function<...>
                access<fext> {};

  using attributes = types::attributes<
      attribute<"mass_matrix", some::matrix<some::scalar, dof, dof>>, q,
      velocity, fext>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) mass_matrix() { return attr<"mass_matrix">(*self()); }

    decltype(auto) velocity() { return Handle::type::velocity::at(*self()); }

    decltype(auto) q() { return Handle::type::q::at(*self()); }

    decltype(auto) fext() { return Handle::type::fext::at(*self()); }
  };
};
}  // namespace siconos::model
