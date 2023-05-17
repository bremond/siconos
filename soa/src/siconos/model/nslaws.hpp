#pragma once

#include "siconos/siconos.hpp"
#include "siconos/utils/pattern.hpp"

namespace siconos {
struct equality_condition_nsl {};
struct relay_nsl {};
struct nonsmooth_law {
  struct newton_impact_friction : item<> {
    struct e : some::scalar, access<e> {};
    struct mu : some::scalar, access<mu> {};

    static constexpr auto size = 2;
    using attributes = types::attributes<e, mu>;

    template <typename Handle>
    struct interface : default_interface<Handle> {
      using default_interface<Handle>::self;
      decltype(auto) e() { return Handle ::type ::e ::at(*self()); }
      decltype(auto) mu() { return Handle ::type ::mu ::at(*self()); }
    };
  };

  struct newton_impact : item<> {
    struct e : some::scalar, access<e>, text<"e"> {};
    static constexpr auto size = 1;
    using attributes = gather<e>;

    template <typename Handle>
    struct interface : default_interface<Handle> {
      using default_interface<Handle>::self;
      decltype(auto) e() { return Handle ::type ::e ::at(*self()); };
    };
  };
};

}  // namespace siconos
