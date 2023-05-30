#pragma once

#include "siconos/utils/pattern.hpp"

namespace siconos::model {

struct lagrangian_r : item<> {
  using attributes = types::attributes<>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

  };
};

  struct lagrangian_tir : item<> {

  };
}  // namespace siconos::model
