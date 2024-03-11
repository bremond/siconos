#pragma once

#include "siconos/simul/simul_head.hpp"

namespace siconos::simul {

template <typename SpaceFilter>
struct interaction_manager : item<> {
  using space_filter = SpaceFilter;
  using nslaw = typename space_filter::nslaw;

  using ncgroups = some::indice_parameter<"ncgroups">;

  using attributes =
      gather<attribute<"space_filter", some::item_ref<space_filter>>,
             attribute<"nslaws", some::matrix<some::item_ref<nslaw>, ncgroups,
                                              ncgroups>>>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) insert_nonsmooth_law(auto nslaw, auto gid1, auto gid2)
    {
      attr<"nslaws"> (*self())(gid1, gid2) = nslaw;
      return *self();
    }

    decltype(auto) get_nonsmooth_law(auto gid1, auto gid2)
    {
      auto& data = self()->data();
      return storage::handle(data, attr<"nslaws">(*self())(gid1, gid2));
    }

    void update_interactions()
    {
      auto& data = self()->data();
      storage::handle(data, attr<"space_filter">(*self())).update_index_set0();
    }

    auto methods()
    {
      using env_t = decltype(self()->env());
      using indice = typename env_t::indice;
      using ref_nslaw_t =
          typename std::decay_t<decltype(attr<"nslaws">(*self()))>::Scalar;

      return collect(
          method("insert_nonsmooth_law",
                 &interface<Handle>::insert_nonsmooth_law<ref_nslaw_t, indice,
                                                          indice>),
          method("get_nonsmooth_law",
                 &interface<Handle>::get_nonsmooth_law<indice, indice>),
          method("update_interactions",
                 &interface<Handle>::update_interactions));
    }
  };

};
}  // namespace siconos::simul
