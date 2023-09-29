#pragma once

#include "siconos/model/model.hpp"

namespace siconos::model {

struct time_invariant {};

namespace lagrangian {
// free functions for all lagrangian relations
template <match::item Relation>
decltype(auto) compute_jachq(Relation& relation)
{
  if constexpr (std::derived_from<Relation, time_invariant>) {
    return relation.h_matrix();
  }
}
}  // namespace lagrangian
template <match::item Nslaw>
struct lagrangian_r : item<> {
  using nslaw = Nslaw;

  using attributes = types::attributes<>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;
  };
};

template <match::item Nslaw>
struct lagrangian_tir
    : item<description<"A lagrangian time invariant relation">>,
      time_invariant {
  using nslaw = Nslaw;
  using nslaw_size = some::indice_value<nslaw::size>;
  using dof = some::indice_parameter<"dof">;

  struct h_matrix : some::matrix<some::scalar, nslaw_size, dof>,
                    access<h_matrix> {};

  struct b : some::scalar, access<b> {};

  using attributes = gather<b, h_matrix>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) h_matrix() { return Handle::type::h_matrix::at(*self()); }

    decltype(auto) b() { return Handle::type::b::at(*self()); }
  };
};

}  // namespace siconos::model
