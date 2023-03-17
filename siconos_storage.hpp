#ifndef SICONOS_STORAGE_HPP
#define SICONOS_STORAGE_HPP

#include <range/v3/detail/variant.hpp>
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
  static constexpr auto memory =
    []<typename T>(typename T::size_type step, T& mem)
    constexpr -> decltype(auto)
  {
    return mem[step%std::size(mem)];
  };

  namespace some
  {
    struct keep : property {};
    struct diagonal : property {};
  }

  struct info{};

  template<match::property ...Parts>
  struct with_properties : item<>
  {
    using properties = gather<Parts...>;
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

  template<match::item T, typename R>
  struct index
  {
    using handle_t = void;
    using type = T;
    using value_t = R;

    value_t _value = {};

    value_t get() const noexcept { return _value; };

    explicit index(R ref) : _value{ref} {}

    index() : _value{} {};

    friend auto operator<=>(const index<T, R>&, const index<T, R>&) = default;

  };

  template<match::item T, typename R, typename D>
  struct handle : index<T, R>, T::template interface<handle<T, R, D>>
  {
    using full_handle_t = void;
    using info_t = std::decay_t<decltype(ground::get<info>(D{}))>;
    using indice = typename info_t::env::indice;
    using attached_storages_t = decltype(filter<hold<decltype(
                                           []<typename X>(X)
      {
        return match::attached_storage<X, T>;
      })>>(typename info_t::all_properties_t{}));
    decltype(auto) attributes() { return typename index<T, R>::type::attributes{}; };

    D& _data;

    template<typename A>
    constexpr decltype(auto) property(A, indice step=0)
    {
      using item_t = T;
      constexpr auto tpl = filter<hold<decltype(
        []<typename X>(X)
        {
          return (match::attached_storage<X, item_t> && (match::tag<X, A>));
        })>>(typename info_t::all_properties_t{});

      static_assert (std::tuple_size_v<decltype(tpl)> >= 1, "attached storage not found"); \

      using attached_storage_t = std::decay_t<decltype(std::get<0>(tpl))>;
      return memory(step, ground::get<attached_storage_t>(data()))[this->get()];
    }

    template<string_literal S>
    constexpr decltype(auto) property( indice step=0)
    {
      return property(symbol<S>{}, step);
    }
    decltype(auto) data() { return _data; };

    explicit handle(R& ref, D& data) : index<T, R>{ref}, _data{data} {};

    explicit handle(index<T, R>& ha, D& data) : index<T, R>{ha}, _data{data} {};

    handle() : index<T, R>{}, _data{} {};

    friend auto operator<=>(const handle<T, R, D>&, const handle<T, R, D>&) = default;

  };


  template<match::item T>
  static decltype(auto) make_half_handle(auto indx)
  {
    using indice = std::decay_t<decltype(indx)>;
    return index<T, indice>{indx};
  }

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

  template<match::property K>
  static auto all_properties_as = [](auto& data) constexpr -> auto
  {
    using info_t = std::decay_t<decltype(ground::get<info>(data))>;
    using all_properties_t = typename info_t::all_properties_t;

    return filter<hold<decltype([]<typename T>(T)
      { return std::derived_from<T, K>; })>>(all_properties_t{});
  };

  template<match::attribute Attr>
  static auto attribute_properties = [](auto& data) constexpr -> auto
  {
    using info_t = std::decay_t<decltype(ground::get<info>(data))>;
    using all_properties_t = typename info_t::all_properties_t;

    return filter<hold<decltype([]<typename T>(T)
      { return std::derived_from<typename T::type, Attr>; })>>(all_properties_t{});
  };

  template<match::attribute Attr, match::property K>
  static auto has_property = [](auto& data) constexpr -> bool
  {
    return
    ground::any_of(
      all_properties_as<K>(data),
      []<match::property P>(P)
      { return std::derived_from<Attr,
                                 typename P::type>; });
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

  template<typename Env, match::attribute ...Attributes>
  struct unit_storage
  {
    using env = Env;

    using type = ground::map<
      ground::pair<Attributes,
                   typename traits::config<env>::template convert<Attributes>::type>
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
      using all_items_t      = decltype(flatten(append(all_items(Items{})...)));
      using all_attributes_t = decltype(flatten(append(all_attributes(Items{})...)));
      using all_properties_t = decltype(flatten(append(all_properties(Items{})...)));

      // subset of all_properties_t
      using all_attached_storages_t = decltype(filter<hold<decltype(
                                                 []<typename T>(T)
        {
          return std::derived_from<T, some::attached_storage>;
        })>>(all_properties_t{}));
    };


    // base map with declared attributes
    using map_t = decltype(
      std::apply(
        []<typename ...Attributes>(Attributes...)
        {
          return typename unit_storage<typename iinfo::env,
                                       Attributes...>::type{};
        },
        flatten(append(typename iinfo::all_attributes_t{},
                       typename iinfo::all_attached_storages_t{}))));

    // final map with info added
    using type = with_info_t<iinfo, map_t>;
  };

  static auto attached_storages = []<typename Handle>(Handle h, auto& data)
    constexpr -> decltype(auto)
  {
    using info_t = std::decay_t<decltype(ground::get<info>(data))>;
    using item_t = typename Handle::type;
    using attached_storages_t = decltype(filter<hold<decltype(
                                          []<typename T>(T)
      {
        return match::attached_storage<T, item_t>;
      })>>(typename info_t::all_properties_t{}));
    return attached_storages_t{};
  };

  template<typename Handle, typename Data>
  using attached_storages_t = std::decay_t<decltype(attached_storages(Handle{}, Data{}))>;

  template<typename T>
  auto make_full_handle(const auto& indx, auto& data)
  {
    using data_t = std::decay_t<decltype(data)>;
    using info_t = std::decay_t<decltype(ground::get<info>(data))>;
    using indice = typename info_t::env::indice;

    indice index = indx;
    return handle<T, indice, data_t>{index, data};
  }

  template<match::item T>
  decltype(auto) make_handle(auto& h)
  {
    return make_full_handle<T>(h.get(), h.data());
  }

  template<typename A>
  static auto get = ground::overload (
    []<match::handle_attribute<A> Handle, typename Data>(
      auto step, Handle& handle, Data& data)
    constexpr -> decltype(auto)
    {
      return memory(step,
                    ground::get<A>(data))
        [handle.get()];
    },
    []<match::handle_attribute<A> Handle, typename Data>(
      Handle& handle, Data& data)
    constexpr -> decltype(auto)
    {
      return memory(0, ground::get<A>(data))
        [handle.get()];
    },
    []<match::handle_attached_storage<A> Handle, typename Data>(
      auto step, Handle& handle, Data& data) constexpr -> decltype(auto)
    {
      using item_t = typename Handle::type;
      using info_t = std::decay_t<decltype(ground::get<info>(data))>;
      constexpr auto tpl = filter<hold<decltype(
        []<typename T>(T)
        {
          return (match::attached_storage<T, item_t> && match::tag<T, A>);
        })>>(typename info_t::all_properties_t{});

//      static_assert (std::tuple_size_v<decltype(tpl)> >= 1, "attached storage not found");
      using attached_storage_t = std::decay_t<decltype(std::get<0>(tpl))>;
      return memory(step, ground::get<attached_storage_t>(data))[handle.get()];
    },
    []<match::handle_attached_storage<A> Handle, typename Data>(
      Handle& handle, Data& data) constexpr -> decltype(auto)
    {
      using item_t = typename Handle::type;
      using info_t = std::decay_t<decltype(ground::get<info>(data))>;
      constexpr auto tpl = filter<hold<decltype(
        []<typename T>(T)
        {
          return (match::attached_storage<T, item_t> && match::tag<T, A>);
        })>>(typename info_t::all_properties_t{});

//      static_assert (std::tuple_size_v<decltype(tpl)> >= 1, "attached storage not found");
      using attached_storage_t = std::decay_t<decltype(std::get<0>(tpl))>;
      return memory(0, ground::get<attached_storage_t>(data))[handle.get()];
    });

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

  static constexpr auto attribute_storage_transform =
    []<typename D, typename F>(D&& d, F&& f)
    constexpr -> decltype(auto)
  {
    return
      ground::map_value_transform(
        std::forward<D>(d),
        [&f]<typename K, typename S>(K, S&& s)
      {
        if constexpr (match::attribute<typename K::type>)
        {
          return std::forward<F>(f)(typename K::type{}, std::forward<S>(s));
        }
        else
        {
          return std::move(s);
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
            (Attribute, Storage& s) -> decltype(auto)
            {
              // if attribute is derived from one of diagonal specifications
              if constexpr (has_property<Attribute, some::diagonal>(base_storage))
              {
                // refine attribute specification toward diagonal storage
                return typename traits::config<typename info_t::env>::template
                  convert<some::diagonal_matrix<typename Attribute::type, Attribute>>::type();
              }
              // if constexpr (has_property<Attribute, other>(base_storage))
              // etc.
              else
              {
                return s;
              }
            }),
          // item level: collection depends on item property
          []<match::item Item, match::attribute Attr>(Item item,
                                                      Attr attr,
                                                      auto s)
          {
            using storage_t = std::decay_t<decltype(s)>;

            if constexpr(std::derived_from<Item,
                         some::undefined_bounded_collection>)
            {
              return ground::pair<Attr,
                typename
                                  Env::template bounded_collection<storage_t, std::get<0>(Item::sizes)>>{}; // std::forward<storage_t>(s)};
            }
            else if constexpr(std::derived_from<Item,
                              some::undefined_unbounded_collection>)
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
        []<match::attribute Attribute>(Attribute attr, auto s)
        {
          using storage_t = std::decay_t<decltype(s)>;
          using all_keeps_t = decltype(all_properties_as<some::keep>(base_storage));
          return memory_t<storage_t, memory_size<Attribute, all_keeps_t>>
            {};
        });
  };

  static auto remove =
    [](auto h, auto& data)
    {
      using item_t = typename std::decay_t<decltype(h)>::type;
      using info_t = std::decay_t<decltype(ground::get<info>(data))>;
      using all_keeps_t = decltype(all_properties_as<some::keep>(data));

      using indice = typename info_t::env::indice;

      auto attrs = flatten(append(attributes(item_t{}),
                                  attached_storages(h, data)));

      using attrs_t = std::decay_t<decltype(attrs)>;

      if constexpr (std::tuple_size_v<attrs_t> > 0)
      {
        ground::for_each(attrs,
                         [&data, &h]<match::attribute A>(A)
                         {
                           return ground::for_each(
                             ground::range<memory_size<A, all_keeps_t>>,
                             [&data, &h](indice step)
                             {
                               move_back(h.get(),
                                         memory(step, ground::get<A>(data)));
                             });
                         });
      }
    };
  template<match::item Item>
  static auto add =
    [](auto&& data) constexpr -> decltype(auto)
    {
      using data_t = std::decay_t<decltype(data)>;
      using info_t = std::decay_t<decltype(ground::get<info>(data))>;
      using all_keeps_t = decltype(all_properties_as<some::keep>(data));

      using indice = typename info_t::env::indice;
      using attached_storage_t = decltype(filter<hold<decltype(
                                           []<typename T>(T)
      {
        return match::attached_storage<T, Item>;
      })>>(typename info_t::all_properties_t{}));

      constexpr auto attrs = flatten(append(attributes(Item{}),
                                            attached_storage_t{}));

      using attrs_t = std::decay_t<decltype(attrs)>;

      // and what about using items = ... ?

      // attributes
      if constexpr (std::tuple_size_v<attrs_t> > 0)
      {
        indice&& index =
          ground::fold_left(attrs, indice{0},
                            [&data]<match::attribute A>(indice k, A)
                            {
                              return ground::fold_left(
                                ground::range<memory_size<A, all_keeps_t>>,
                                k,
                                [&data](indice n, auto step)
                                {
                                  auto&& storage = memory(step, ground::get<A>(std::forward<data_t>(data)));
                                  using storage_t = std::decay_t<decltype(storage)>;
                                  if constexpr (match::push_back<storage_t>)
                                  {
                                    storage.push_back(typename storage_t::value_type{});
                                    return std::size(storage) - 1;
                                  }
                                  else
                                  {
                                    storage[0] = typename storage_t::value_type{};
                                    return 0;
                                  }
                                });
                            });
        return make_full_handle<Item>(index, data);
      }
      else
      {
        return make_full_handle<Item>(indice{0}, data);
      }
    };

  template<match::item T>
  static auto for_each_attribute = [](auto&& fun, auto& data)
    constexpr
  {
    ground::for_each(_attributes(T{}),
      [&fun]<match::attribute ...Attrs>(Attrs&...)
      {
        (fun(Attrs{}), ...);
      });

  };

  template<typename T>
  struct access
  {
    static constexpr auto at =
      ground::overload(
        []<
        typename Data,
        typename U = T,
        match::handle<decltype(item_attribute<U>(typename std::decay_t<decltype(ground::get<info>(Data{}))>::all_items_t{}))> Handle
        >
        (Handle h, Data& data)
        -> decltype(auto)
        {
          return siconos::get<U>(h, data);
        },
        []<
        typename U = T,
        typename FullHandle
        >
        (FullHandle h)
        -> decltype(auto)
        {
          return siconos::get<U>(h, h.data());
        });
      };

  static auto for_each = [](auto m, auto&& fun) constexpr -> void
  {
    ground::for_each(
      m,
      ground::dup(ground::lockstep(fun)(ground::first, ground::second)));
  };
}
#endif
