#pragma once

#include "siconos/utils/environment.hpp"
#include "siconos/utils/some.hpp"

namespace siconos {
struct aaa {
  int v = 1;
};
struct bbb : item<> {
  struct attr : some::specific<pointer<aaa>>,
                siconos::access<attr>,
                some::attribute<> {};
  using attributes = gather<attr>;
  template <typename H>
  struct interface : siconos::default_interface<H> {
  };
};

static_assert(
    std::is_same_v<
        traits::config<standard_environment>::convert<some::scalar>::type,
        double>);
static_assert(
    std::is_same_v<
        traits::config<standard_environment>::convert<pointer<float>>::type,
        pointer<float>>);

static_assert(std::is_same_v<traits::config<standard_environment>::convert<
                                 some::specific<pointer<float>>>::type,
                             pointer<float>>);

}  // namespace siconos
