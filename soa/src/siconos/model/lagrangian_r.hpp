#pragma once

#include <functional>  // mem_fn

#include "siconos/algebra/algebra.hpp"
#include "siconos/model/lagrangian_ds.hpp"
#include "siconos/model/model.hpp"
namespace siconos::storage::pattern::match {
template <typename T>
concept linear_relation = match::handle<T, model::linear>;

template <typename T>
concept relation1 = match::handle<T, model::relation1>;

template <typename T>
concept relation2 = match::handle<T, model::relation2>;

}  // namespace siconos::storage::pattern::match

namespace siconos::model {

template <auto NSLSize>
struct lagrangian_r : item<>,
                      linear,
                      relation1,
                      relation2,
                      any_lagrangian_relation {
  using nslaw_size = some::param_val<NSLSize>;
  using dof = some::indice_parameter<"dof">;

  using attributes = gather<
      attribute<"b", some::scalar>,
      attribute<"h_matrix", some::matrix<some::scalar, nslaw_size, dof>>>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) h_matrix() { return attr<"h_matrix">(*self()); }

    decltype(auto) b() { return attr<"b">(*self()); }

    decltype(auto) compute_jachq(auto step, auto& q1, auto& q2,
                                 auto& h_matrix1, auto& h_matrix2)
    {
      h_matrix1 = h_matrix();
      h_matrix2 = -h_matrix();
    }
    decltype(auto) compute_jachq(auto step, auto& q, auto& h_matrix1)
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
