#pragma once

#include "siconos/storage/ground/ground.hpp"

namespace siconos::storage {

// the storage info key
struct info {};

// pre-map cases
template <typename... Pairs>
auto get_info(ground::tuple<Pairs...>&& data)
{
  return ground::second(ground::at(
      static_cast<ground::tuple<Pairs...>&&>(data), ground::size_c<0>));
}

template <typename... Pairs>
auto get_info(const ground::tuple<Pairs...>& data)
{
  return ground::second(ground::at(data, ground::size_c<0>));
}

// template <typename... Pairs>
// auto get_info(ground::tuple<Pairs...>& data)
// {
//   return ground::second(ground::at(data, ground::size_c<0>));
// }

// database case
template <typename... Pairs>
auto get_info(ground::database<Pairs...>&& data)
{
  return ground::get<info>(static_cast<ground::database<Pairs...>&&>(data));
}

template <typename D>
using get_info_t = decltype(get_info(std::decay_t<D>{}));

}  // namespace siconos::storage
