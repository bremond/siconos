#pragma once

#include "siconos/siconos.hpp"
#include "siconos/utils/pattern.hpp"

namespace siconos {
struct lcp {};

template <typename Type>
struct one_step_nonsmooth_problem : item<> {
  using problem_type = Type;
  struct level : some::indice {};
  using attributes = gather<level>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) level() { return Handle ::type ::level ::at(*self()); };
  };
};
}  // namespace siconos
