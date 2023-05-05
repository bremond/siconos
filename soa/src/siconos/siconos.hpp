#pragma once

#include "siconos/storage/storage.hpp"
#include "siconos/utils/pattern.hpp"
#include "siconos/algebra/linear_algebra.hpp"
#include "siconos/algebra/numerics.hpp"
#include <range/v3/view/zip.hpp>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/filter.hpp>
#include <tuple>

#define GET(X) decltype(auto) X()                                       \
  { return Handle::type::X::at(*self()); }

namespace siconos
{

  decltype(auto) edge1(auto& g, auto& descr)
  {
    auto [oei, oee] = g.out_edges(descr);
    return oei;
  };

  decltype(auto) edge2(auto& g, auto& descr)
  {
    auto [oei, oee] = g.out_edges(descr);
    return oee;
  };



}

