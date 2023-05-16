#pragma once

#include "siconos/siconos.hpp"
#include "siconos/utils/pattern.hpp"

namespace siconos {
template <typename Nslaw, typename Relation, std::size_t K = 2>
struct interaction : item<> {
  using dof = some::indice_parameter<"dof">;
  using nslaw_size = some::indice_value<Nslaw::size>;

  struct nonsmooth_law : some::item_ref<Nslaw>, access<nonsmooth_law> {};
  struct relation : some::item_ref<Relation>, access<relation> {};

  struct h_matrix : some::matrix<some::scalar, nslaw_size, dof>,
                    access<h_matrix> {};
  struct lambda : some::vector<some::scalar, nslaw_size>, access<lambda> {};

  struct y : some::vector<some::scalar, nslaw_size>, access<y> {};

  using attributes = gather<dof, nslaw_size, nonsmooth_law, relation, h_matrix, lambda, y>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) nonsmooth_law()
    {
      return Handle ::type ::nonsmooth_law ::at(*self());
    }
    decltype(auto) lambda() { return Handle ::type ::lambda ::at(*self()); }
    decltype(auto) relation()
    {
      return handle(Handle::type::relation::at(*self()), self()->data());
    }

    decltype(auto) h_matrix() { return Handle::type::h_matrix::at(*self()); }

    decltype(auto) y() { return Handle ::type ::y ::at(*self()); }
  };
};

}  // namespace siconos
