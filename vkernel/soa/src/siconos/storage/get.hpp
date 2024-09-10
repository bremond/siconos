#pragma once

#include "siconos/storage/ground/ground.hpp"
#include "siconos/storage/pattern/base.hpp"
#include "siconos/storage/pattern/base_concepts.hpp"
#include "siconos/storage/pattern/pattern.hpp"

namespace siconos::storage {

using namespace pattern;

template <typename A>
static auto get = ground::overload(
    // get<Attr>(data, step, handle)
    []<match::handle_attribute<A> Handle, typename Data>(
        Data&& data, auto step, Handle&& handle) constexpr -> decltype(auto) {
      return memory(step,
                    ground::get<A>(static_cast<Data&&>(data)))[handle.get()];
    },
    // get<Attr>(data, step, handle)
    []<match::handle_attribute<A> Handle, typename Data>(
        Data& data, auto step, Handle& handle) constexpr -> decltype(auto) {
      return memory(step, ground::get<A>(data))[handle.get()];
    },
    // get<Attr>(data, handle)
    []<match::handle_attribute<A> Handle, typename Data>(
        Data& data, Handle& handle) constexpr -> decltype(auto) {
      return memory(0, ground::get<A>(data))[handle.get()];
    }  // ,
       // // get<Attached_storage>(data, step, h)
       // []<match::handle_attached_storage<A> Handle, typename Data>(
       //     Data& data, auto step, Handle& handle) constexpr ->
       //     decltype(auto)
       //     {
       //   using item_t = typename Handle::type;
       //   using info_t = std::decay_t<decltype(ground::get<info>(data))>;
       //   constexpr auto tpl = ground::filter(
       //       typename info_t::all_properties_t{},
       //       ground::is_a_model<[]<typename T>() consteval {
       //         return (match::attached_storage<T, item_t> &&
       //         match::tag<T, A>);
       //       }>);

    //   //      static_assert (std::tuple_size_v<decltype(tpl)> >= 1,
    //   "attached
    //   //      storage not found");
    //   using attached_storage_t =
    //   std::decay_t<decltype(std::get<0>(tpl))>; return memory(step,
    //                 ground::get<attached_storage_t>(data))[handle.get()];
    // },
    // // get<Attached_storage>(data,h)
    // []<match::handle_attached_storage<A> Handle, typename Data>(
    //   Data& data, Handle& handle) constexpr -> decltype(auto) {
    //   using item_t = typename Handle::type;
    //   using info_t = std::decay_t<decltype(ground::get<info>(data))>;

    //   constexpr auto tpl = ground::filter(
    //       typename info_t::all_properties_t{},
    //       ground::is_a_model<[]<typename T>() consteval {
    //         return (match::attached_storage<T, item_t> &&
    //         match::tag<T, A>);
    //       }>);
    //   //      constexpr auto tpl = filter<hold<decltype([]<typename
    //   T>(T) {
    //   //        return (match::attached_storage<T, item_t> &&
    //   match::tag<T,
    //   //        A>);
    //   //      })>>(typename info_t::all_properties_t{});

    //   //      static_assert (std::tuple_size_v<decltype(tpl)> >= 1,
    //   "attached
    //   //      storage not found");
    //   using attached_storage_t =
    //   std::decay_t<decltype(std::get<0>(tpl))>; return memory(0,
    //   ground::get<attached_storage_t>(data))[handle.get()];
    // }
);

template <typename T>
struct access {
  static constexpr auto at = ground::overload(
      []<typename Data, typename U = T,
         match::handle<decltype(item_attribute<U>(
             typename get_info_t<Data>::all_items_t{}))>
             Handle>(Handle h, Data& data) -> decltype(auto) {
        return siconos::storage::get<U>(h, data);
      },
      []<typename U = T, typename FullHandle>(FullHandle h)
          -> decltype(auto) { return siconos::storage::get<U>(h.data(), h); },
      []<typename U = T, typename FullHandle>(
          FullHandle h, typename FullHandle::indice step) -> decltype(auto) {
        return siconos::storage::get<U>(h.data(), step, h);
      });
};

template <string_literal S>
static auto param = [](auto h) constexpr -> decltype(auto) {
  return h.template env_param<S>();
};

template <string_literal S>
static auto prop = [](auto h) constexpr -> decltype(auto) {
  return h.template property<S>();
};

template <string_literal S>
static auto attr = []<typename H>(H h, typename H::indice step =
                                           0) constexpr -> decltype(auto) {
  using attr_n = attr_t<typename H::type, S>;
  return memory(step, ground::get<attr_n>(h.data()))[h.get()];
};

template <match::attribute T>
static constexpr decltype(auto) attr_memory(auto& data)
{
  return ground::get<T>(data);
};

template <match::item I, string_literal S>
static constexpr decltype(auto) attr_memory(auto& data)
{
  return ground::get<attr_t<I, S>>(data);
};

template <string_literal S>
static constexpr auto is_identified_by =
    ground::is_a_model<[]<typename T>() constexpr {
      return match::tag<T, symbol<S>>;
    }>;

template <match::item I, string_literal S>
static auto prop_memory = [](auto& data) constexpr -> decltype(auto) {
  using info_t = get_info_t<decltype(data)>;
  constexpr auto tpl =
      ground::filter(ground::filter(typename info_t::all_properties_t{},
                                    is_attached_storage<I>),
                     is_identified_by<S>);

  //  constexpr auto tpl = filter<hold<decltype([]<typename X>(X) {
  //      return (match::attached_storage<X, item_t> && (match::tag<X,
  //      symbol<S>>));
  //    })>>(typename info_t::all_properties_t{});

  static_assert(ground::size(tpl) >= ground::size_c<1_c>,
                "attached storage not found");

  using attached_storage_t = std::decay_t<decltype(tpl[0_c])>;
  return ground::get<attached_storage_t>(data);
};

template <match::attribute T>
static constexpr decltype(auto) attr_values(auto& data, auto step)
{
  return memory(step, (attr_memory<T>(data)));
};

template <match::item I, string_literal S, typename D>
static constexpr decltype(auto) attr_values(D&& data, auto step)
{
  return memory(step, (attr_memory<I, S>(data)));
};

template <match::item I, string_literal S>
static auto prop_values =
    [](auto& data, auto step) constexpr -> decltype(auto) {
  return memory(step, (prop_memory<I, S>(data)));
};

// fix H is a handle defined in storage
template <typename Hc>
static auto constexpr methods(Hc hc)
{
  using handle_t = typename decltype(+hc)::type;
  auto data = typename handle_t::data_t{};
  auto h = add<typename handle_t::type>(data);

  if constexpr (match::methods<handle_t>) {
    return h.methods();
  }
  else {
    return gather<>{};
  }
}

template <typename Astor>
static auto constexpr attached_storage_name(Astor astor)
{
  using tag = typename Astor::tag;
  if constexpr (match::symbol<tag>) {
    return tag{};
  }
  else {
    []<bool flag = false>() {
      ground::type_trace<tag>();
      static_assert(flag, "tag is not a symbol!");
    }();
  }
}
}  // namespace siconos::storage
