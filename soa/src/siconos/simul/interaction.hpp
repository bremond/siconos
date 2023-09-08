#pragma once

#include "siconos/simul/simul.hpp"

namespace siconos::simul {

template <typename Nslaw, typename Relation, std::size_t K = 2>
struct interaction : item<> {
  using nslaw_t = Nslaw;
  using relation_t = Relation;

  using dof = some::indice_parameter<"dof">;
  using nslaw_size = some::indice_value<Nslaw::size>;

  struct nonsmooth_law : some::item_ref<Nslaw>, access<nonsmooth_law> {};
  struct relation : some::item_ref<Relation>, access<relation> {};

  struct h_matrix1 : some::matrix<some::scalar, nslaw_size, dof>,
                     access<h_matrix1> {};

  struct h_matrix2 : some::matrix<some::scalar, nslaw_size, dof>,
                     access<h_matrix2> {};

  struct lambda : some::vector<some::scalar, nslaw_size>, access<lambda> {};
  struct y : some::vector<some::scalar, nslaw_size>, access<y> {};
  struct ydot : some::vector<some::scalar, nslaw_size>, access<ydot> {};

  using attributes = gather<dof, nslaw_size, nonsmooth_law, relation,
                            h_matrix1, h_matrix2, lambda, y, ydot>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) nonsmooth_law()
    {
      return Handle ::type ::nonsmooth_law ::at(*self());
    }

    decltype(auto) relation()
    {
      return make_handle<relation_t>(*self());
    }

    decltype(auto) h_matrix1() { return Handle::type::h_matrix1::at(*self()); }
    decltype(auto) h_matrix2() { return Handle::type::h_matrix2::at(*self()); }

    decltype(auto) lambda() { return Handle ::type ::lambda ::at(*self()); }
    decltype(auto) y() { return Handle ::type ::y ::at(*self()); }
    decltype(auto) ydot() { return Handle ::type ::y ::at(*self()); }
  };
};

}  // namespace siconos::simul
