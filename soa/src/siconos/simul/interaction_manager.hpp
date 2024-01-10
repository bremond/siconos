#pragma once

#include "siconos/simul/simul.hpp"

namespace siconos::simul {

template <typename Nslaw>
struct interaction_manager : item<> {
  using ncgroups = some::indice_parameter<"ncgroups">;

  using attributes = gather<attribute<
      "nslaws", some::matrix<some::item_ref<Nslaw>, ncgroups, ncgroups>>>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) insert_nonsmooth_law(auto nslaw, auto gid1, auto gid2)
    {
      attr<"nslaws">(*self())(gid1, gid2) = nslaw;
      return *self();
    }

    void update_interactions()
    {
      // broad phase
    }

    auto methods()
    {
      using env_t = decltype(self()->env());
      using indice = typename env_t::indice;
      using ref_nslaw_t = typename std::decay_t<decltype(attr<"nslaws">(*self()))>::Scalar;

      return collect(method("insert_nonsmooth_law",
                            &interface<Handle>::insert_nonsmooth_law<ref_nslaw_t, indice ,indice>),
                     method("update_interactions",
                            &interface<Handle>::update_interactions));
    }
  };
};
}  // namespace siconos::simul
