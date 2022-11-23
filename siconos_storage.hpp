#ifndef SICONOS_STORAGE_HPP
#define SICONOS_STORAGE_HPP

#include <type_traits>
#include <variant>
#include <concepts>
#include <cassert>

#include "siconos_ground.hpp"
#include "siconos_util.hpp"
#include "siconos_traits.hpp"
#include "siconos_pattern.hpp"

namespace siconos
{
  namespace some
  {
    struct keep : kind {};
    struct diagonal : kind {};
  }

  struct info{};

  template<match::kind ...Parts>
  struct with_kinds : item<>
  {
    using kinds = gather<Parts...>;
  };

  template<match::attribute Attr, std::size_t N>
  struct keep : some::keep
  {
    using type = Attr;
    using keep_t = void;
//    using attribute = Attr;
    static constexpr std::size_t size = N;
  };

  template<match::matrix M>
  struct diagonal : some::diagonal
  {
    using type = M;
    using diagonal_t = void;
  };

  template<match::item I>
  struct unbounded_collection : I, some::unbounded_collection
  {
    using type = I;
  };

  template<match::item I, std::size_t N>
  struct bounded_collection : I, some::bounded_collection
  {
    using type = I;
    static constexpr std::size_t size = N;
  };

  template<typename T, std::size_t N>
  using memory_t = std::array<T, N>;

  template<match::attribute Attr, concepts::tuple_like keeps_t>
  static constexpr std::size_t memory_size = []()
  {
    // rewrite this with a filter...
    // let rec loop(...) = ...
    constexpr auto rec_loop = []<typename Tpl>(auto&& loop, Tpl) constexpr
    {
      constexpr auto keep = car(Tpl{});
      using keep_t = std::decay_t<decltype(keep)>;

      if constexpr (std::is_same_v<Attr, typename keep_t::type>)
      {
        return keep_t::size;
      }
      else if constexpr (std::tuple_size_v<Tpl> > 1)
      {
        return loop(loop, cdr(Tpl{}));
      }
      else
      {
        // memory size not specified
        return 1;
      }
    };

    if constexpr (std::tuple_size_v<keeps_t> > 0)
    {
      return rec_loop(rec_loop, keeps_t{});
    }
    else
    {
      // no keep in memory has been specified
      return 1;
    }
  } ();

  static constexpr auto all_kinds =
    rec([]<match::item Item>(auto&& all_kinds, Item)
  {
    auto ckinds = []() constexpr
    {
      if constexpr (match::kinds<Item>)
      {
        return typename Item::kinds{};
      }
      else
      {
        return std::tuple<>{};
      }
    }();

    if constexpr (match::items<Item>)
    {
      return flatten(append(
                       ckinds,
                       transform(all_kinds,
                                 typename Item::items{})));
    }
    else
    {
      return ckinds;
    }
  });

  template<match::kind K>
  static auto all_kinds_of = [](auto& data) constexpr -> auto
  {
    using info_t = std::decay_t<decltype(ground::get<info>(data))>;
    using all_kinds_t = typename info_t::all_kinds_t;

    return filter<hold<decltype([]<typename T>(T)
      { return std::derived_from<T, K>; })>>(all_kinds_t{});
  };

  template<match::attribute Attr>
  static auto attribute_kinds = [](auto& data) constexpr -> auto
  {
    using info_t = std::decay_t<decltype(ground::get<info>(data))>;
    using all_kinds_t = typename info_t::all_kinds_t;

    return filter<hold<decltype([]<typename T>(T)
      { return std::derived_from<typename T::type, Attr>; })>>(all_kinds_t{});
  };

  template<match::attribute Attr, match::kind K>
  static auto is_kind_of = [](auto& data) constexpr -> bool
  {
    return
    ground::any_of(
      all_kinds_of<K>(data),
      []<match::kind D>(D)
      { return std::derived_from<Attr,
                                 typename D::type>; });
  };

  static constexpr auto memory =
    []<typename T>(typename T::size_type step, T& mem)
    constexpr -> decltype(auto)
  {
    return mem[step%std::size(mem)];
  };

  template<match::attribute T>
  static auto get_memory = [](auto& data) constexpr
    -> decltype(auto)
  {
    return ground::get<T>(data);
  };

  static auto move_back = [](const auto i, auto& a) constexpr
  {
    if constexpr (match::push_back<std::decay_t<decltype(a)>>)
    {
      a[i] = std::move(a.back());
      a.pop_back();
    }
    // else...
  };

  template<typename Env, typename ...Attributes>
  struct unit_storage
  {
    using env = Env;

    using type = ground::map<
      ground::pair<Attributes,
                   typename traits::config<env, Attributes>::type>
      ...>;

  };

  template<typename Info, typename M>
  using with_info_t = decltype(ground::insert(M{}, ground::pair<info, Info>{}));

  template<typename Env, match::item ...Items>
  struct item_storage
  {
    struct iinfo
    {
      using env = Env;
      using all_items_t = decltype(flatten(append(all_items(Items{})...)));
      using all_attributes_t = decltype(flatten(append(all_attributes(Items{})...)));
      using all_kinds_t = decltype(flatten(append(all_kinds(Items{})...)));
    };

    using map_t = decltype(
      std::apply(
        []<typename ...Attributes>(Attributes...)
        {
          return typename unit_storage<typename iinfo::env,
                                       Attributes...>::type{};
        },
        typename iinfo::all_attributes_t{}));

    using type = with_info_t<iinfo, map_t>;
  };

  template<match::attribute A>
  static auto get = ground::overload (
    []<match::handle_attribute<A> Handle, typename Data>(
      auto step, Handle handle, Data&& data)
    constexpr -> decltype(auto)
    {
      return memory(step,
                    ground::get<A>(std::forward<std::decay_t<Data>>(data)))
        [handle.get()];
    },
    []<match::handle_attribute<A> Handle, typename Data>(
      Handle handle, Data&& data)
    constexpr -> decltype(auto)
    {
      return memory(0, ground::get<A>(std::forward<std::decay_t<Data>>(data)))
        [handle.get()];
    },
    []<typename Data>(
      Data&& data)
    constexpr -> decltype(auto)
    {
      return ground::get<A>(std::forward<std::decay_t<Data>>(data))[0];
    }
    );


  static constexpr auto item_storage_transform =
    []<typename D>(D&& d, auto&&f) constexpr -> decltype(auto)
  {
    using info_t = std::decay_t<decltype(ground::get<info>(d))>;

    return
    ground::map_transform(
      d,
      [&f]<typename P>(P&& key_value)
      {
        static auto&& key = ground::first(key_value);
        static auto&& value = ground::second(key_value);
        using key_t = std::decay_t<decltype(key)>;
        using value_t = std::decay_t<decltype(value)>;

        if constexpr (match::attribute<typename key_t::type>)
        {
          using attr_t = typename key_t::type;
          return f(item_attribute<attr_t>(typename info_t::all_items_t{}),
                   attr_t{},
                   std::forward<value_t>(value));
        }
        else
        {
          return std::forward<std::decay_t<P>>(key_value);
        }
      });
  };

  static constexpr auto attribute_storage_transform = [](auto&& d, auto&& f)
    constexpr -> decltype(auto)
  {
    return
      ground::map_value_transform(
      d,
      [&f]<typename K, typename S>(K, S&& s)
      {
        if constexpr (match::attribute<typename K::type>)
        {
          return f(typename K::type{}, s);
        }
        else
        {
          return std::forward<S>(s);
        }
      });
  };

  template<typename Env, match::item ...Items>
  static auto make_storage = []()
    constexpr -> decltype(auto)
  {
    auto base_storage = typename item_storage<Env, Items...>::type{};
    using info_t = std::decay_t<decltype(ground::get<info>(base_storage))>;
    return
      attribute_storage_transform(
        item_storage_transform(
          attribute_storage_transform(
            base_storage,
            // attribute level for base storage specifications
            [&base_storage]<match::attribute Attribute, typename Storage>
            (Attribute, Storage&& s)
            {
              // if attribute is derived from one of diagonal specifications
              if constexpr (is_kind_of<Attribute, some::diagonal>(base_storage))
              {
                // refine attribute specification toward diagonal storage
                return typename traits::config<typename info_t::env,
                                               some::diagonal_matrix<Attribute>>::type();
              }
              else
              {
                return std::forward<std::decay_t<Storage>>(s);
              }
            }),
          // item level: collection depends on item kind
          []<match::item Item, match::attribute Attr>(Item item,
                                                      Attr attr,
                                                      auto&& s)
          {
            using storage_t = std::decay_t<decltype(s)>;

            if constexpr(std::derived_from<Item, some::bounded_collection>)
            {
              return ground::pair<Attr,
                typename
                Env::template bounded_collection<storage_t, Item::size>>{}; // std::forward<storage_t>(s)};
            }
            else if constexpr(std::derived_from<Item, some::unbounded_collection>)
            {
              return ground::pair<Attr,
                typename Env:: template unbounded_collection<storage_t>>{};
            }
            else // default storage
            {
              return ground::pair<Attr, typename Env:: template default_storage<storage_t>>{};
            }
          }),

        // attribute level: memory depends on keeps
        []<match::attribute Attribute>(Attribute attr, auto&& s)
        {
          using storage_t = std::decay_t<decltype(s)>;
          using all_keeps_t = decltype(all_kinds_of<some::keep>(base_storage));
          return memory_t<storage_t, memory_size<Attribute, all_keeps_t>>
            {};
        });
  };

  template<match::item Item>
  static auto add =
    [](auto&& data) constexpr -> decltype(auto)
    {
      using data_t = std::decay_t<decltype(data)>;
      using info_t = std::decay_t<decltype(ground::get<info>(data))>;
      using all_keeps_t = decltype(all_kinds_of<some::keep>(data));

      using indice = typename info_t::env::indice;

      constexpr auto attrs = attributes(Item{});
      using attrs_t = std::decay_t<decltype(attrs)>;

      // and what about using items = ... ?

      // attributes
      if constexpr (std::tuple_size_v<attrs_t> > 0)
      {
        return make_internal_handle<Item>(
          ground::fold_left(attributes(Item{}), indice{0},
                            [&data]<match::attribute A>(indice k, A)
                            {
                              return ground::fold_left(
                                ground::range<memory_size<A, all_keeps_t>>,
                                k,
                                [&data](indice n, auto step)
                                {
                                  auto& storage = memory(step, ground::get<A>(std::forward<data_t>(data)));
                                  using storage_t = std::decay_t<decltype(storage)>;
                                  if constexpr (match::push_back<storage_t>)
                                  {
                                    storage.push_back(typename storage_t::value_type{});
                                    return n + std::size(storage) - 1;
                                  }
                                  else
                                  {
                                    // same place
                                    storage[step] = typename storage_t::value_type{};
                                    return n;
                                  }
                                });
                            }) / std::tuple_size_v<attrs_t>
          );
      }
      else
      {
        return make_internal_handle<Item>(empty{});
      }
    };

  template<match::item T>
  static auto for_each_attribute = [](auto&& fun, auto& data)
    constexpr
  {
    ground::for_each(all_attributes(T{}),
      [&fun]<match::attribute ...Attrs>(Attrs&...)
      {
        (fun(Attrs{}), ...);
      });

  };

  template<match::item T>
  static auto add_attributes = [](auto& data) constexpr -> decltype(auto)
  {
    using data_t = std::decay_t<decltype(data)>;
    using env = typename data_t::env;
    using indice = typename env::indice;

    return
    ground::fold_left(all_attributes(T{}), indice{0},
                      [&data]<match::attribute A>(indice n, A)
                      {
                        auto step = get_step<A>(data);
                        auto& mem = get_memory<A>(data);
                        auto& storage = memory(step, mem);

                        using storage_t = std::decay_t<decltype(storage)>;
                        if constexpr (match::push_back<storage_t>)
                        {
                          // append new place
                          storage.push_back(typename storage_t::value_type{});
                          return n + std::size(storage) - 1;
                        }
                        else
                        {
                          // same place
                          storage[0] = typename storage_t::value_type{};
                          return n;
                        }
                      }) / std::tuple_size_v<decltype(all_attributes(T{}))>;
  };

  template<typename T>
  struct access
  {
    static constexpr auto internal_get =
      []<typename U = T>(match::internal_handle auto h, auto& data)
      -> decltype(auto)
    {
      return siconos::get<U>(h, data);
    };

    static constexpr auto at = []<typename U = T>(match::handle auto h)
      -> decltype(auto)
    {
      return fix_map(siconos::get<U>)(h);
    };

  };

  template<typename T>
  static constexpr auto fixed_add = fix(add<T>);

}
#endif
