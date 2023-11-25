#pragma once

#include <cassert>
#include <concepts>
#include <range/v3/detail/variant.hpp>
#include <tuple>
#include <type_traits>
#include <variant>

#include "siconos/storage/ground/ground.hpp"
#include "siconos/storage/pattern/base.hpp"
#include "siconos/storage/pattern/base_concepts.hpp"
#include "siconos/storage/pattern/pattern.hpp"
#include "siconos/storage/some/some.hpp"
#include "siconos/storage/traits/traits.hpp"
#include "siconos/utils/range.hpp"

namespace siconos::storage {

using namespace pattern;

template <match::item I, typename Tag, match::attribute DataSpec>
struct attached : DataSpec, some::attached_storage {
  using item = I;
  using tag = Tag;
  using data_spec = DataSpec;
};

static constexpr auto memory = []<typename T>(
                                   typename std::decay_t<T>::size_type step,
                                   T&& mem) constexpr -> decltype(auto) {
  return static_cast<T&&>(mem)[step % std::size(static_cast<T&&>(mem))];
};

namespace property {
struct keep : some::property {};

struct time_invariant : some::property {};

struct refine : some::property {};

struct diagonal : refine {
  template <match::attribute A>
  using refine = some::diagonal_matrix<typename A::type, A>;
};
struct unbounded_diagonal : refine {
  template <match::attribute A>
  using refine = some::unbounded_diagonal_matrix<typename A::type>;
};
}  // namespace property

template <match::attribute A, typename T>
struct refine_with_type : A {
  using type = T;
};

struct info {};

template <match::property... Parts>
struct with_properties : item<> {
  using properties = gather<Parts...>;
};

template <match::attribute Attr, std::size_t N>
struct keep : property::keep {
  using type = Attr;
  using keep_t = void;
  //    using attribute = Attr;
  static constexpr std::size_t size = N;
};

template <match::attribute Attr>
struct time_invariant : property::time_invariant {
  using type = Attr;
  using time_invariant_t = void;
};

//  template <typename T, string_literal = "">
// struct diagonal {
//};

template <match::item Item, string_literal S>
struct diagonal : property::diagonal {
  using type = attr_t<Item, S>;
  static_assert(match::abstract_matrix<type>);
  using diagonal_t = void;
};

// template <match::abstract_matrix M>
// struct diagonal<M> : property::diagonal {
//   using type = M;
//   using diagonal_t = void;
// };

template <match::abstract_matrix M>
struct unbounded_diagonal : property::unbounded_diagonal {
  using type = M;
  using unbounded_diagonal_t = void;
};

template <match::item T, typename R>
struct index {
  using handle_t = void;
  using type = T;
  using value_t = R;

  value_t _value = {};

  value_t get() const noexcept { return _value; };

  explicit index(R ref) : _value{ref} {}

  index() : _value{} {};

  friend auto operator<=>(const index<T, R>&, const index<T, R>&) = default;
};

template <typename Handle>
struct default_interface {
  decltype(auto) self()
  {
    return static_cast<Handle*>(
        this);  // handle inherits from default_interface
  };
};

template <match::item T, typename R, typename D>
struct handle : index<T, R>, T::template interface<handle<T, R, D>> {
  using full_handle_t = void;
  using info_t = std::decay_t<decltype(ground::get<info>(D{}))>;
  using indice = typename info_t::env::indice;
  using attached_storages_t =
      decltype(filter<hold<decltype([]<typename X>(X) {
        return match::attached_storage<X, T>;
      })>>(typename info_t::all_properties_t{}));

  D& _data;

  constexpr decltype(auto) attributes()
  {
    return attributes(typename index<T, R>::type{});
  };

  template <typename A>
  constexpr decltype(auto) property(A, indice step = 0)
  {
    using item_t = T;
    constexpr auto tpl = filter<hold<decltype([]<typename X>(X) {
      return (match::attached_storage<X, item_t> && (match::tag<X, A>));
    })>>(typename info_t::all_properties_t{});

    static_assert(std::tuple_size_v<decltype(tpl)> >= 1,
                  "attached storage not found");

    using attached_storage_t = std::decay_t<decltype(std::get<0>(tpl))>;
    return memory(step, ground::get<attached_storage_t>(data()))[this->get()];
  }

  // not convenient, it needs to specify template keyword:
  // some_handle.template property<S>() = ...
  template <string_literal S>
  constexpr decltype(auto) property(indice step = 0)
  {
    return property(symbol<S>{}, step);
  }

  decltype(auto) data() { return _data; };

  explicit handle(D& data, R& ref) : index<T, R>{ref}, _data{data} {};

  explicit handle(D& data, index<T, R>& ha) : index<T, R>{ha}, _data{data} {};

  explicit handle(D& data, index<T, R>&& ha)
      : index<T, R>{ha}, _data{data} {};

  handle() : index<T, R>{}, _data{} {};

  handle(const handle& h) : index<T, R>(h.get()), _data(h._data){};
  handle(handle&&) = default;
  handle operator=(const handle& h) { return handle(h); };
  handle operator=(handle&& h) { return handle(h); };
  ;
  friend auto operator<=>(const handle<T, R, D>&,
                          const handle<T, R, D>&) = default;
};

template <match::item T>
static decltype(auto) make_half_handle(auto indx)
{
  using indice = std::decay_t<decltype(indx)>;
  return index<T, indice>{indx};
}

template <typename T, std::size_t N>
using memory_t = std::array<T, N>;

template <match::attribute Attr, concepts::tuple_like keeps_t>
static constexpr std::size_t memory_size = []() {
  // rewrite this with a filter...
  // let rec loop(...) = ...
  constexpr auto rec_loop = []<typename Tpl>(auto&& loop, Tpl) constexpr {
    constexpr auto keep = car(Tpl{});
    using keep_t = std::decay_t<decltype(keep)>;

    if constexpr (std::is_same_v<Attr, typename keep_t::type>) {
      return keep_t::size;
    }
    else if constexpr (std::tuple_size_v<Tpl> > 1) {
      return loop(loop, cdr(Tpl{}));
    }
    else {
      // memory size not specified
      return 1;
    }
  };

  if constexpr (std::tuple_size_v < keeps_t >> 0) {
    return rec_loop(rec_loop, keeps_t{});
  }
  else {
    // no keep in memory has been specified
    return 1;
  }
}();

template <match::property K>
static auto all_properties_as = [](auto& data) constexpr -> auto
{
  using info_t = std::decay_t<decltype(ground::get<info>(data))>;
  using all_properties_t = typename info_t::all_properties_t;

  return ground::filter(all_properties_t{}, ground::derive_from<K>);
};

template <match::attribute Attr>
static auto attribute_properties = [](auto& data) constexpr -> auto
{
  using info_t = std::decay_t<decltype(ground::get<info>(data))>;
  using all_properties_t = typename info_t::all_properties_t;

  return filter<hold<decltype([]<typename T>(T) {
    return std::derived_from<typename T::type, Attr>;
  })>>(all_properties_t{});
};

template <match::attribute Attr, match::property K>
static constexpr bool has_property(auto& data)
{
  return ground::any_of(all_properties_as<K>(data), []<match::property P>(P) {
    return std::derived_from<Attr, typename P::type>;
  });
}

template <match::item Item, match::property K>
static constexpr auto has_property(auto& data)
{
  return ground::any_of(typename Item::properties{},
                        []<match::property P>(P) {
                          return std::derived_from<typename P::type, K>;
                        });
}

template <typename A, typename K, typename D>
using has_property_t = std::decay_t<decltype(has_property<A, K>(D{}))>;

static auto refine_attribute = []<match::attribute Attr, typename D>(
                                   const D& data,
                                   Attr) constexpr -> decltype(auto) {
  using refines = decltype(filter<hold<decltype([]<match::property P>(P) {
    return std::derived_from<Attr, typename P::type>;
  })>>(all_properties_as<property::refine>(data)));
  if constexpr (std::tuple_size_v < refines >> 0) {
    return typename nth_t<0, refines>::template refine<Attr>{};
  }
  else {
    // return attribute as it is
    return Attr{};
  }
};

static constexpr auto refine_recursively_attribute =
    []<match::attribute Attr>(auto& data, Attr) {
      using data_t = std::decay_t<decltype(data)>;
      constexpr auto rec_loop = []<typename IAttr>(auto&& loop,
                                                   IAttr) constexpr {
        if constexpr (match::attribute_with_internal_type<IAttr>) {
          // look inside and apply fun
          if constexpr (match::attribute<typename IAttr::type>) {
            using r =
                refine_with_type<IAttr, decltype(loop(
                                            loop, typename IAttr::type{}))>;
            return refine_attribute(data_t{}, r{});
          }
          else {
            return refine_attribute(data_t{}, IAttr{});
          }
        }
        else {
          // terminal attribute
          return refine_attribute(data_t{}, IAttr{});
        }
      };

      return rec_loop(rec_loop, Attr{});
    };

static auto move_back = [](const auto i, auto& a) constexpr {
  if constexpr (match::push_back<std::decay_t<decltype(a)>>) {
    a[i] = std::move(a.back());
    a.pop_back();
  }
  // else...
};

template <typename Env, match::attribute... Attributes>
struct unit_storage {
  using env = Env;

  using type = ground::map<ground::pair<
      Attributes,
      typename traits::config<env>::template convert<Attributes>::type>...>;
};

template <typename Info, typename M>
using with_info_t = decltype(ground::insert(M{}, ground::pair<info, Info>{}));

template <typename Env, match::item... Items>
struct item_storage {
  struct iinfo {
    using env = Env;
    using all_items_t = decltype(ground::std_tuple(ground::fold_left(
        // take last defined items
        ground::reverse(flatten(append(all_items(Items{})...))),
        ground::make_tuple(), [](auto acc, auto elem) {
          if constexpr (ground::contains(acc, elem)) {
            return acc;
          }
          else {
            return ground::append(acc, elem);
          };
        })));
    using all_attributes_t =
        decltype(flatten(append(all_attributes(Items{})...)));
    using all_properties_t =
        decltype(flatten(append(all_properties(Items{})...)));

    // subset of all_properties_t
    using all_attached_storages_t =
        decltype(filter<hold<decltype([]<typename T>(T) {
          return std::derived_from<T, some::attached_storage>;
        })>>(all_properties_t{}));
  };

  // base map with declared attributes
  using map_t = decltype(std::apply(
      []<typename... Attributes>(Attributes...) {
        return
            typename unit_storage<typename iinfo::env, Attributes...>::type{};
      },
      flatten(append(typename iinfo::all_attributes_t{},
                     typename iinfo::all_attached_storages_t{}))));

  // final map with info added
  using type = with_info_t<iinfo, map_t>;
};

static auto attached_storages =
    []<typename Handle>(Handle h, auto& data) constexpr -> decltype(auto) {
  using info_t = std::decay_t<decltype(ground::get<info>(data))>;
  using item_t = typename Handle::type;
  using attached_storages_t =
      decltype(filter<hold<decltype([]<typename T>(T) {
        return match::attached_storage<T, item_t>;
      })>>(typename info_t::all_properties_t{}));
  return attached_storages_t{};
};

template <typename Handle, typename Data>
using attached_storages_t =
    std::decay_t<decltype(attached_storages(Handle{}, Data{}))>;

template <typename T>
auto make_full_handle(auto& data, const auto& indx)
{
  using data_t = std::decay_t<decltype(data)>;
  using info_t = std::decay_t<decltype(ground::get<info>(data))>;
  using indice = typename info_t::env::indice;

  indice index = indx;
  return handle<T, indice, data_t>{data, index};
}

template <match::item T>
decltype(auto) make_handle(auto& h)
{
  return make_full_handle<T>(h.data(), h.get());
}

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
    //     Data& data, auto step, Handle& handle) constexpr -> decltype(auto)
    //     {
    //   using item_t = typename Handle::type;
    //   using info_t = std::decay_t<decltype(ground::get<info>(data))>;
    //   constexpr auto tpl = ground::filter(
    //       typename info_t::all_properties_t{},
    //       ground::is_a_model<[]<typename T>() consteval {
    //         return (match::attached_storage<T, item_t> && match::tag<T,
    //         A>);
    //       }>);

    //   //      static_assert (std::tuple_size_v<decltype(tpl)> >= 1,
    //   "attached
    //   //      storage not found");
    //   using attached_storage_t = std::decay_t<decltype(std::get<0>(tpl))>;
    //   return memory(step,
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
    //         return (match::attached_storage<T, item_t> && match::tag<T,
    //         A>);
    //       }>);
    //   //      constexpr auto tpl = filter<hold<decltype([]<typename T>(T) {
    //   //        return (match::attached_storage<T, item_t> && match::tag<T,
    //   //        A>);
    //   //      })>>(typename info_t::all_properties_t{});

    //   //      static_assert (std::tuple_size_v<decltype(tpl)> >= 1,
    //   "attached
    //   //      storage not found");
    //   using attached_storage_t = std::decay_t<decltype(std::get<0>(tpl))>;
    //   return memory(0,
    //   ground::get<attached_storage_t>(data))[handle.get()];
    // }
);

static constexpr auto item_storage_transform =
    []<typename D>(D&& d, auto&& f) constexpr -> decltype(auto) {
  using info_t = std::decay_t<decltype(ground::get<info>(d))>;

  return ground::map_transform(d, [&f]<typename P>(P&& key_value) {
    static auto&& key = ground::first(key_value);
    static auto&& value = ground::second(key_value);
    using key_t = std::decay_t<decltype(key)>;
    using value_t = std::decay_t<decltype(value)>;

    if constexpr (match::attribute<typename key_t::type>) {
      using attr_t = typename key_t::type;
      return f(item_attribute<attr_t>(typename info_t::all_items_t{}),
               attr_t{}, std::forward<value_t>(value));
    }
    else {
      return std::forward<std::decay_t<P>>(key_value);
    }
  });
};

static constexpr auto attribute_storage_transform =
    []<typename D, typename F>(D&& d, F&& f) constexpr -> decltype(auto) {
  return ground::map_value_transform(
      std::forward<D>(d), [&f]<typename K, typename S>(K, S&& s) {
        if constexpr (match::attribute<typename K::type>) {
          return std::forward<F>(f)(typename K::type{}, std::forward<S>(s));
        }
        else {
          return std::move(s);
        }
      });
};

template <typename Env, match::item... Items>
static auto make = []() constexpr -> decltype(auto) {
  auto base_storage = typename item_storage<Env, Items...>::type{};
  using info_t = std::decay_t<decltype(ground::get<info>(base_storage))>;
  return attribute_storage_transform(
      item_storage_transform(
          attribute_storage_transform(
              base_storage,
              // attribute level for base storage specifications
              []<match::attribute Attribute, typename Storage>(
                  Attribute, Storage& s) -> decltype(auto) {
                // if attribute is derived from one of diagonal
                // specifications, etc.
                return typename traits::config<typename info_t::env>::
                    template convert<decltype(refine_recursively_attribute(
                        base_storage, Attribute{}))>::type{};
              }),
          // item level: collection depends on item property
          []<match::item Item, match::attribute Attr>(Item item, Attr attr,
                                                      auto s) {
            using storage_t = std::decay_t<decltype(s)>;

            if constexpr (match::wrap<Item>) {
              return ground::pair<
                  Attr,
                  typename traits::config<Env>::template convert<
                      typename Item::template wrapper<storage_t>>::type>{};
            }
            else  // default storage
            {
              return ground::pair<
                  Attr, typename Env::template default_storage<storage_t>>{};
            }
          }),
      // attribute level: memory depends on keeps
      []<match::attribute Attribute>(Attribute attr, auto s) {
        using storage_t = std::decay_t<decltype(s)>;
        using all_keeps_t =
            decltype(all_properties_as<property::keep>(base_storage));
        return memory_t<storage_t, memory_size<Attribute, all_keeps_t>>{};
      });
};

static auto remove = [](auto& data, auto& h) {
  using item_t = typename std::decay_t<decltype(h)>::type;
  using info_t = std::decay_t<decltype(ground::get<info>(data))>;
  using all_keeps_t = decltype(all_properties_as<property::keep>(data));

  using indice = typename info_t::env::indice;

  auto attrs =
      flatten(append(attributes(item_t{}), attached_storages(h, data)));

  using attrs_t = std::decay_t<decltype(attrs)>;

  if constexpr (std::tuple_size_v < attrs_t >> 0) {
    ground::for_each(attrs, [&data, &h]<match::attribute A>(A) {
      return ground::for_each(ground::range<memory_size<A, all_keeps_t>>,
                              [&data, &h](indice step) {
                                move_back(h.get(),
                                          memory(step, ground::get<A>(data)));
                              });
    });
  }
};

template <match::item Item>
static auto add = [](auto&& data) constexpr -> decltype(auto) {
  using data_t = std::decay_t<decltype(data)>;
  using info_t = std::decay_t<decltype(ground::get<info>(data))>;
  using all_keeps_t = decltype(all_properties_as<property::keep>(data));

  using indice = typename info_t::env::indice;
  auto attached_storage =
      ground::filter(typename info_t::all_properties_t{},
                     ground::is_a_model<[]<typename T>() consteval {
                       return match::attached_storage<T, Item>;
                     }>);

  constexpr auto attrs =
      flatten(append(attributes(Item{}), attached_storage));

  using attrs_t = std::decay_t<decltype(attrs)>;

  // and what about using items = ... ?

  // attributes
  if constexpr (std::tuple_size_v < attrs_t >> 0) {
    indice&& index = ground::fold_left(
        attrs, indice{0}, [&data]<match::attribute A>(indice k, A) {
          return ground::fold_left(
              ground::range<memory_size<A, all_keeps_t>>, k,
              [&data](indice n, auto step) {
                auto&& storage =
                    memory(step, ground::get<A>(std::forward<data_t>(data)));
                using storage_t = std::decay_t<decltype(storage)>;
                if constexpr (match::push_back<storage_t>) {
                  storage.push_back(typename storage_t::value_type{});
                  return std::size(storage) - 1;
                }
                else {
                  storage[0] = typename storage_t::value_type{};
                  return 0;
                }
              });
        });
    return make_full_handle<Item>(data, index);
  }
  else {
    return make_full_handle<Item>(data, indice{0});
  }
};

template <match::item T>
static constexpr void for_each_attribute(T)
{
  return ground::compose(ground::for_each, attributes)(T{});
};

template <typename T>
struct access {
  static constexpr auto at = ground::overload(
      []<typename Data, typename U = T,
         match::handle<decltype(item_attribute<U>(
             typename std::decay_t<decltype(ground::get<info>(
                 Data{}))>::all_items_t{}))>
             Handle>(Handle h, Data& data)
          ->decltype(auto) { return siconos::storage::get<U>(h, data); },
      []<typename U = T, typename FullHandle>(FullHandle h)->decltype(auto) {
        return siconos::storage::get<U>(h.data(), h);
      },
      []<typename U = T, typename FullHandle>(
          FullHandle h, typename FullHandle::indice step)
          ->decltype(auto) {
            return siconos::storage::get<U>(h.data(), step, h);
          });
};

// ?
// static auto for_each = [](auto m, auto&& fun) constexpr -> void {
//  ground::for_each(
//      m, ground::dup(ground::lockstep(fun)(ground::first, ground::second)));
//};

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

template <match::item I, string_literal S>
static auto prop_memory = [](auto& data) constexpr -> decltype(auto) {
  using info_t = std::decay_t<decltype(ground::get<storage::info>(data))>;
  using item_t = I;
  constexpr auto tpl =
      ground::filter(typename info_t::all_properties_t{},
                     ground::is_a_model<[]<typename T>() consteval {
                       return (match::attached_storage<T, item_t> &&
                               match::tag<T, symbol<S>>);
                     }>);

  //  constexpr auto tpl = filter<hold<decltype([]<typename X>(X) {
  //      return (match::attached_storage<X, item_t> && (match::tag<X,
  //      symbol<S>>));
  //    })>>(typename info_t::all_properties_t{});

  static_assert(std::tuple_size_v<decltype(tpl)> >= 1,
                "attached storage not found");

  using attached_storage_t = std::decay_t<decltype(std::get<0>(tpl))>;
  return ground::get<attached_storage_t>(data);
};

template <match::attribute T>
static constexpr decltype(auto) attr_values(auto& data, auto step)
{
  return memory(step, (attr_memory<T>(data)));
};

template <match::item I, string_literal S>
static constexpr decltype(auto) attr_values(auto& data, auto step)
{
  return memory(step, (attr_memory<I, S>(data)));
};

template <match::item I, string_literal S>
static auto prop_values =
    [](auto& data, auto step) constexpr -> decltype(auto) {
  return memory(step, (prop_memory<I, S>(data)));
};

template <match::item I>
static auto handles = [](auto& data, auto step) constexpr -> decltype(auto) {
  using info_t = std::decay_t<decltype(ground::get<storage::info>(data))>;
  using env = typename info_t::env;
  using indice = typename env::indice;

  using attributes_t = std::decay_t<decltype(attributes(I{}))>;
  // need at least one attributes
  // empty items are not supposed to exist
  indice num = std::size(attr_values<nth_t<0, attributes_t>>(data, step));
  return views::iota((indice)0, num) | views::transform([&data](indice i) {
           return handle<I, indice, std::decay_t<decltype(data)>>(
               data, index<I, indice>(i));
         });
};

/*siconos::storage::ground::dump_keys(data, [](auto&& s) { std::cout << s<<
 * std::endl;});*/

using pattern::attr_t;
using pattern::wrap;

}  // namespace siconos::storage

namespace siconos {
using siconos::storage::access;
using siconos::storage::default_interface;
using siconos::storage::prop;
}  // namespace siconos
