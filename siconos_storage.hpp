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
  template<typename T, std::size_t N>
  using memory_t = std::array<T, N>;

  template<match::attribute Attr, concepts::tuple_like keeps_t>
  static constexpr std::size_t memory_size = []()
  {
    // let rec loop(...) = ...
    constexpr auto rec_loop = []<typename Tpl>(auto&& loop, Tpl) constexpr
    {
      constexpr auto keep = car(Tpl{});
      using keep_t = std::decay_t<decltype(keep)>;

      if constexpr (std::is_same_v<Attr, typename keep_t::attribute>)
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

    return rec_loop(rec_loop, keeps_t{});
  } ();

  static constexpr auto all_keeps =
    rec([]<match::item Item>(auto&& all_keeps, Item)
  {
    auto ckeeps = []() constexpr
    {
      if constexpr (match::keeps<Item>)
      {
        return typename Item::keeps{};
      }
      else
      {
        return std::tuple<>{};
      }
    }();

    if constexpr (match::items<Item>)
    {
      return flatten(append(
                       ckeeps,
                       transform(all_keeps,
                                 typename Item::items{})));
    }
    else
    {
      return ckeeps;
    }
  });

  static constexpr auto memory =
    []<typename T>(typename T::size_type step, T& mem)
    constexpr -> typename T::value_type&
  {
    return mem[step%std::size(mem)];
  };

  template<match::attribute T>
  static auto get_memory = [](auto& data) constexpr
    -> decltype(auto)
  {
    using data_t = std::decay_t<decltype(data)>;
    constexpr typename data_t::data_types t = T{};
    auto& mem = std::get<t.index()>(data._collections);
    return mem;
  };

  template<match::attribute T>
  static auto get_step = [](auto& data) constexpr -> decltype(auto)
  {
    using data_t = std::decay_t<decltype(data)>;
    using time_discretization = typename data_t::machine::time_discretization;

    if constexpr
      (memory_size<T, typename data_t::machine::osi::keeps> > 1)
    {
      // one time_discretization / data
      return get_memory<typename time_discretization::step>(data)[0][0];
    }
    else
    {
      return 0;
    }
  };

  template<match::attribute T>
  static auto get =
    [](match::internal_handle auto ihandle, auto& data)
    constexpr -> decltype(auto)
    {
      using data_t = std::decay_t<decltype(data)>;
      using topology_t = typename data_t::topology_t;
      using item_t = decltype(item_attribute<T>(typename data_t::all_items_t{}));
      using indice = typename data_t::indice;

      auto& mem = get_memory<T>(data);
      indice step = get_step<T>(data);

      if constexpr (match::vertex_item<item_t>)
      {
        auto& dsg = get_memory<typename topology_t::dsg0>(data)[0][0];
        auto indx = dsg.index(ihandle.get());
        return memory(step, mem)[indx];
      }
      else
      {
        return memory(step, mem)[ihandle.get()];
      }
    };

  template<typename Env, match::attribute T, typename Items>
  static auto attribute_storage =
      []() constexpr
      {
        using item_t = decltype(item_attribute<T>(Items{}));

        static_assert(match::item<item_t>);

        if constexpr (match::graph_item<item_t>)
        {
          // a collection
          return
            (typename Env::template collection<
             typename traits::config<Env, T>::type>{});
        }
        else
        {
          static_assert(!match::vertex_item<item_t> && !match::edge_item<item_t>);
          return (typename place_holder<typename traits::config<Env, T>::type>::type{});
        }
      };

  template<typename Env, match::attribute T, typename Items>
  using attribute_storage_t = decltype(attribute_storage<Env, T, Items>());


/*  template<match::attribute T>
  static auto set =
    []<typename H>(H handle)
    constexpr -> decltype(auto)
    {
      auto& data = handle.data();
      using data_t = std::decay_t<decltype(data)>;

      return overload {
        [handle, &data](typename attribute_storage_t<typename data_t::env, T,
                        typename data_t::all_items_t>::value_type&& value)
      {
        using value_t = std::decay_t<decltype(value)>;
        using item_t = decltype(item_attribute<T>(typename data_t::all_items_t{}));
        using indice = typename data_t::indice;

        auto& mem = get_memory<T>(data);
        indice step = get_step<T>(data);

        if constexpr (match::vertex_item<item_t>)
        {
          auto indx = data._graph.index(handle.get());
          return memory(step, mem)[indx] = std::forward<value_t>(value);
        }
        else
        {
          return memory(step, mem)[handle.get()] = std::forward<value_t>(value);
        }
      },
      [handle, &data](typename attribute_storage_t<typename data_t::env, T,
                      typename data_t::all_items_t>::value_type& value)
      {
        using item_t = decltype(item_attribute<T>(typename data_t::all_items_t{}));
        using indice = typename data_t::indice;

        auto& mem = get_memory<T>(data);
        indice step = get_step<T>(data);

        if constexpr (match::vertex_item<item_t>)
        {
          auto indx = data._graph.index(handle.get());
          return memory(step, mem)[indx] = std::move(value);
        }
        else
        {
          return memory(step, mem)[handle.get()] = std::move(value);
        }
      },
      []<bool flag = false, typename X>(X) { static_assert(flag, "set: unknown type");}
      };
    };
*/
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

  template<typename Env, typename ...Items>
  struct item_storage
  {
    using env = Env;
    using all_items_t = decltype(flatten(all_items(Items{})...));
    using attributes = decltype(flatten(all_attributes(Items{})...));

    using items_t = decltype(std::apply([]<typename ...IItems>(IItems...)
      {
        return ground::map<
          ground::pair<IItems,
                       decltype(all_attributes(IItems{}))>
          ...>{};
      }, all_items_t{}));

    using storage_t = decltype(std::apply([]<typename ...Attributes>(Attributes...)
                                          { return typename unit_storage<env, Attributes...>::type{} ;}, attributes{}));

    items_t xitems;
    storage_t storage;
  };

  static auto iiget = ground::overload
    (
      []<match::attribute A, typename D>(A, D& data) constexpr -> decltype(auto)
      { return ground::get<A>(data.storage); },

      []<match::item I, typename D>(I, D& data) constexpr -> decltype(auto)
      { return ground::get<I>(data.xitems); });

  template<typename T>
  static auto iget = ground::t_arg<T>(iiget);

  static constexpr auto storage_transform =
    []<typename D>(D&& d, auto&&f) constexpr -> decltype(auto)
  {
    using data_t = std::decay_t<decltype(d)>;
    return
    ground::map_value_transform(
      d.storage,
      [&f]<typename K, typename S>(K, S&& store)
      {
        using A = typename K::type;
        return f(item_attribute<A>(typename data_t::all_items_t{}), std::forward<S>(store));
      });
  };

  namespace xmemory
  {
    template<typename T, std::size_t N>
    using wrap = std::array<T, N>;

    static constexpr auto get = [](auto&& mem, auto step)
      constexpr -> decltype(auto)
    {
      return mem[step%std::size(mem)];
    };
  }

  template<typename Env, typename Mach, typename ...Attributes>
  struct storage
  {
    using env = Env;
    using machine = Mach;
    using topology_t = typename machine::topology;
    using all_items_t = decltype(all_items(Mach{}));
    using all_keeps_t = decltype(all_keeps(Mach{}));

    using indice = typename env::indice;
    using scalar = typename env::scalar;
    using graph_t = typename env::graph;

    static_assert(std::tuple_size_v<all_items_t> > 0);
    static_assert((match::attribute<Attributes> && ...));

    // wrap collection type in a memory container according to
    // osi::keeps specification
    template<typename A>
    using collection =
      memory_t<attribute_storage_t<env, A, all_items_t>,
               memory_size<A, all_keeps_t>>;

    // global variant for all attributes
    using data_types = std::variant<Attributes ...>;

    // global storage type for all defined collections
    using collections_t = std::tuple<collection<Attributes>...>;

    // global storage
    collections_t _collections;

    // administrative counter
    indice _counter = 0;

    ///////////////////////
    // methods & operators
    template<typename Fun>
    decltype (auto) operator () (Fun&& fun)
    {
      return proj(*this)(fun);
    }

  };

  template<typename Sim>
  struct machinery
  {
    using vertex_items =
      decltype(
        std::apply(
          []<typename ...Vertex_items>(Vertex_items...)
          {
            return std::variant<Vertex_items ...>{};
          }, all_vertex_items(Sim{})));

    struct vertex : vertex_item<>
    {
      using time_discretization = typename Sim::time_discretization;
      using osi = typename Sim::one_step_integrator;
      using topology = typename Sim::topology;
      struct descriptor : some::vdescriptor<vertex_items> {};

      using attributes = gather<descriptor>;
      using items = gather<Sim>;
    };

  };

  template<typename Env, match::item Sim>
  static auto make_storage = []() constexpr
  {
    using machine = machinery<Sim>;

    return std::apply(
      []<typename ...Attributes>(Attributes...)
      {
        return storage<Env, typename machine::vertex, Attributes...>{};
      }, all_attributes(typename machine::vertex{}));
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

  template<match::vertex_item T>
  static auto add_vertex_item = [](auto& data) constexpr -> decltype(auto)
  {
    using data_t = std::decay_t<decltype(data)>;
    using topology_t = typename data_t::topology_t;

    using vdescriptor = typename data_t::machine::vertex::descriptor;
    //    auto step = 0 ; //get_step<T>(data);

    // auto& dsg = get<topology_t::dsg0>(handle?, data);
    // only one topology
    auto& dsg = get_memory<typename topology_t::dsg0>(data)[0][0];

    auto vd = dsg.add_vertex(data._counter++);

    get_memory<vdescriptor>(data)[0].push_back({vd, T{}});

    dsg.index(vd) = std::size(get_memory<vdescriptor>(data)[0])-1;

    [[maybe_unused]] auto idx = add_attributes<T>(data);

    assert(idx == dsg.index(vd));

    return make_internal_handle<T>(vd);
  };

  template<match::item T>
  static auto add = [](auto& data) constexpr -> decltype(auto)
  {
    if constexpr (match::vertex_item<T>)
    {
      return add_vertex_item<T>(data);
    }
    else
    {
      return make_internal_handle<T>(add_attributes<T>(data));
    }
  };

  static auto remove_vertex_item = [](auto handle,
                                      const auto step) constexpr ->decltype(auto)
  {
    auto& data = handle.data();
    using data_t = std::decay_t<decltype(data)>;
    using topology_t = typename data_t::topology_t;
    using vdescriptor = typename data_t::machine::vertex::descriptor;
    auto& vd_store = get_memory<vdescriptor>(data)[0];

    //auto& dsg = get<topology_t::dsg0>(data);
    auto& dsg = get_memory<typename topology_t::dsg0>(data)[0][0];

    auto i = dsg.index(handle.get());

    dsg.remove_vertex(handle.get());

    move_back(i, vd_store);
    dsg.index(std::get<0>(vd_store[i])) = i;

    std::visit([&i,&data]<typename Item>(Item)
    {
      for_each_attribute<Item>(
        [&i,&data]<match::attribute A>(A &a)
        {
          auto step = get_step<A>(data);
          auto& mem = get_memory<A>(data);
          auto& storage = memory(step, mem);

          move_back(i, storage);
        },
        data);
    }, (std::get<1>(vd_store[i])));
    return 0;
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

    static constexpr auto get = []<typename U = T>(match::handle auto h)
      -> decltype(auto)
    {
      return fix_map(siconos::get<U>)(h);
    };

  };

  template<typename T>
  static constexpr auto fixed_add = fix(add<T>);
}
#endif
