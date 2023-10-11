#pragma once

#include "siconos/storage/storage.hpp"
#include "siconos/storage/pattern/pattern.hpp"
#include "siconos/storage/some/some.hpp"

namespace siconos::model
{
  using namespace storage;
  using namespace storage::pattern;

  struct linear : some::property {};
  struct any_lagrangian_relation {};

}
