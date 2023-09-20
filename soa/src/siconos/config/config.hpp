#pragma once

#include "siconos/storage/storage.hpp"

namespace siconos::config
{
  template<typename ...Cfs>
  using map = storage::ground::map<Cfs...>;

  template<typename K, typename V>
  using assoc = storage::ground::pair<K, V>;

  using storage::pattern::param;

  using storage::pattern::param_val;

  template<storage::pattern::string_literal S, auto N>
  using iparam = assoc<param<S>, param_val<N>>;
}
