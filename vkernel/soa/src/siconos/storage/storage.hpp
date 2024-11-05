#pragma once

#include "siconos/utils/range.hpp"
#include <boost/hana/string.hpp>

import siconos.storage;
import siconos.storage.ground;

namespace siconos::storage {

static auto apply_fun = []<typename Item, typename SomeFun>(
                            auto& data, Item, SomeFun&& some_fun) {
  using item_t = Item;
  using info_t = get_info_t<decltype(data)>;
  using all_keeps_t = decltype(all_properties_as<property::keep>(data));

  using indice = typename info_t::env::indice;

  auto attrs = ground::tuple_unique(
      concat(attributes(item_t{}), attached_storages(item_t{}, data)));

  if constexpr (ground::size(attrs) > ground::size_c<0>) {
    ground::for_each(attrs, [&data, &some_fun]<match::attribute A>(A) {
      return ground::for_each(ground::range<memory_size<A, all_keeps_t>()>,
                              [&data, &some_fun](indice step) {
                                static_cast<SomeFun&&>(some_fun)(
                                    memory(step, ground::get<A>(data)));
                              });
    });
  }
};



template <match::item I>
constexpr decltype(auto) handles(auto& data, std::size_t step = 0) {
  using info_t = get_info_t<decltype(data)>;
  using env = typename info_t::env;
  using indice = typename env::indice;

  using attributes_t = std::decay_t<decltype(attributes(I{}))>;
  // need at least one attributes
  // empty items are not supposed to exist
  indice num =
      std::size(attr_values<ground::nth_t<0, attributes_t>>(data, step));
  return view::iota((indice)0, num) | view::transform([&data](indice i) {
           return handle<I, indice, std::decay_t<decltype(data)>>(
               data, index<I, indice>(i));
         });
};


using pattern::attr_t;
using pattern::wrap;

using namespace boost::hana::literals;

/* the neighborhood engine */

}  // namespace siconos::storage
