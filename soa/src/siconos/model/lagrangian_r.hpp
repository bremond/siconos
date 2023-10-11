#pragma once

#include "siconos/model/lagrangian_ds.hpp"
#include "siconos/model/model.hpp"
#include "siconos/storage/pattern/pattern.hpp"


namespace siconos::storage::pattern::match {
template <typename T>
concept linear_relation =
match::handle_property<T, model::linear> &&
  match::handle<T, model::any_lagrangian_relation>;
}


namespace siconos::model {

template <match::item Nslaw>
struct lagrangian_r : item<>, any_lagrangian_relation {
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

    decltype(auto) compute_jachq(auto step, auto& h_matrix1, auto& h_matrix2)
    {
      h_matrix1 = h_matrix();
      h_matrix2 = -h_matrix();
    }
    decltype(auto) compute_jachq(auto step, auto& h_matrix1)
    {
      h_matrix1 << -h_matrix();
    }
  };
};

template <match::item Nslaw>
struct diskdisk_r : item<>, any_lagrangian_relation {
  using nslaw = Nslaw;
  using nslaw_size = some::indice_value<nslaw::size>;
  using dof = some::indice_parameter<"dof">;

  struct h_matrix : some::matrix<some::scalar, nslaw_size, dof>,
                    access<h_matrix> {};

  struct b : some::scalar, access<b> {};

  using attributes = gather<b,
                            attribute<"r1", some::scalar>,
                            attribute<"r2", some::scalar>, h_matrix>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) r1() { return attr<"r1">(*self()); };
    decltype(auto) r2() { return attr<"r2">(*self()); };

    decltype(auto) h_matrix() { return Handle::type::h_matrix::at(*self()); }

    decltype(auto) b() { return Handle::type::b::at(*self()); }

    decltype(auto) compute_jachq(auto step, auto& h_matrix1, auto& h_matrix2, auto& q1, auto&q2)
    {
      h_matrix1 = h_matrix();
      h_matrix2 = -h_matrix();
    }
    decltype(auto) compute_jachq(auto step, auto& h_matrix1)
    {
      h_matrix1 << -h_matrix();
    }
  };
};


namespace lagrangian {
// free functions for all lagrangian relations


template <typename Data, match::linear_relation HandleRel,
          match::handle<lagrangian_ds> HandleDS>
decltype(auto) compute_h(Data& data, HandleRel& relation, HandleDS& ds)
{
  return relation.h_matrix() * ds.q();
}
template <typename Data, match::item Relation>
decltype(auto) compute_jachq(Data& data, Relation& relation)
{
  if constexpr (has_property_t<Relation, property::time_invariant, Data>{}) {
    return relation.h_matrix();
  }
}
}  // namespace lagrangian

}  // namespace siconos::model
