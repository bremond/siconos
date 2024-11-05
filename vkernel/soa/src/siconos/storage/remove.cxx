module;

#include <type_traits>
#include <algorithm>
#include <assert.h>

export module siconos.storage:remove;

import siconos.storage.ground;
import siconos.storage.pattern;

import :info;
import :properties;
import :memory;

export namespace siconos::storage {

using namespace pattern;

auto move_back = [](const auto i, auto& a) constexpr {
  if constexpr (match::push_back<std::decay_t<decltype(a)>>) {
    assert((int)a.size() >= 1);

    a[i] = std::move(a.back());
    a.pop_back();

    assert((int)a.size() >= 0);
  }
  // else...
};

auto remove = [](auto& data, auto& h) {
  using item_t = typename std::decay_t<decltype(h)>::type;
  using info_t = get_info_t<decltype(data)>;
  using all_keeps_t = decltype(all_properties_as<property::keep>(data));

  using indice = typename info_t::env::indice;

  auto attrs = ground::tuple_unique(
      concat(attributes(item_t{}), attached_storages(h.item_type(), data)));

  if constexpr (ground::size(attrs) > ground::size_c<0>) {
    ground::for_each(attrs, [&data, &h]<match::attribute A>(A) {
      return ground::for_each(ground::range<memory_size<A, all_keeps_t>()>,
                              [&data, &h](indice step) {
                                move_back(h.get(),
                                          memory(step, ground::get<A>(data)));
                              });
    });
  }
};

}  // namespace siconos::storage
