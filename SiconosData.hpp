#ifndef SICONOS_DATA_HPP
#define SICONOS_DATA_HPP

#include <type_traits>
#include <variant>
#include <tuple>
#include <fmt/core.h>
#include <fmt/ranges.h>
#include <boost/pfr.hpp>
#include <boost/hana.hpp>
#include <boost/hana/ext/std/tuple.hpp>

using fmt::print;


namespace siconos
{
  struct any
  {
    struct scalar{};
    struct indice{};
    struct vdescriptor{};

    template<size_t N, size_t M>
    struct matrix{};
    template<size_t N>
    struct vector{};
  };


  namespace traits
  {
    template<typename T>
    struct missing_conversion_for{};

    template<typename E, typename T>
    struct config
    {
      using type = missing_conversion_for<T>;
    };

    template<typename E>
    struct config<E, any::scalar>
    {
      using type = typename E::scalar;
    };

    template<typename E>
    struct config<E, any::indice>
    {
      using type = typename E::indice;
    };

    template<typename E>
    struct config<E, any::vdescriptor>
    {
      using type = typename E::vdescriptor;
    };

    template<std::size_t N, typename E>
    struct config<E, any::vector<N>>
    {
      using type = typename E::template vector<N>;
    };

    template<std::size_t N, std::size_t M, typename E>
    struct config<E, any::matrix<N, M>>
    {
      using type = typename E::template matrix<N, M>;
    };
  };

  template<typename U, typename... T>
  constexpr bool contains(std::tuple<T...>)
  {
    return (std::is_same_v<U, T> || ...);
  }

  auto proj = [] (auto& data)
  {
    return
      ([&data](auto&& fun) -> decltype(auto)
      {
        return
          ([&data,&fun]<typename ...As>(As&&...args) -> decltype(auto)
          {
            return fun(std::forward<As>(args)..., data);
          });
      });
  };

  auto compose = [](auto&& f, auto&& g) -> decltype(auto)
  {
    return
      [&f,&g]<typename ...As>(As&& ...args) -> decltype(auto)
      {
        return f(g(std::forward<As>(args)...));
      };
  };

  template<typename T>
  constexpr auto get_array = [](auto& data) -> decltype(auto)
  {
    using data_t = std::decay_t<decltype(data)>;
    constexpr typename data_t::data_types t = T{};
    auto& array = std::get<t.index()>(data._collections);
    return array;
  };

  template<typename T>
  constexpr auto get = [](const auto vd, const auto step, auto& data) -> decltype(auto)
  {
    auto indx = data._graph.index(vd);
    auto& array = get_array<T>(data);
    return array[step % std::size(array)][indx];
  };

  void move_back(const auto i, auto& a)
  {
    a[i] = std::move(a.back());
    a.pop_back();
  }

  template<typename Env, typename OSI, typename Item>
  static constexpr auto make_item_data();

  template<typename Env, typename OSI, typename ...Items>
  struct data
  {
    using osi = OSI;
    using env = Env;
    using indice = typename env::indice;
    using scalar = typename env::scalar;
    using graph_t = typename env::graph;

    template<typename C>
    using collection =
      typename std::conditional<contains<C>(typename osi::keep::type{}),
                                std::array<typename env::template collection<typename traits::config<env, typename C::type>::type>, osi::keep::value>,
                                std::array<typename env::template collection<typename traits::config<env, typename C::type>::type>, 1>>::type;

    template<typename T>
    struct item_access
    {
      template<typename ...Ds>
      struct item_data
      {
        using types = std::variant<Ds...>;
        using collections_t = std::tuple<collection<Ds>...>;

      };

      using attributes = typename T::attributes;

      using collections_t = decltype(std::apply([]<typename ...Ts>(Ts&&...)
                                                {
                                                  return typename item_data<Ts...>::collections_t{};
                                                },
                                                attributes{}));
    };

    using item_types = std::variant<Items...>;
//    std::tuple<item_access<Items>...> _items;

    using tpl_cat = decltype(std::tuple_cat<typename Items::attributes...>
                             (typename Items::attributes{}...));

    using data_types = decltype(std::apply(
                                  []<typename ...Ts>(Ts&&...)
                                  {
                                    return std::variant<Ts...>{};
                                  },
                                  std::declval<tpl_cat>()));

    using collections_t = decltype(std::tuple_cat<typename item_access<Items>::collections_t...>
                                   (typename item_access<Items>::collections_t{}...));

    /* stored data */
    indice _counter = 0;
    graph_t _graph;
    collections_t _collections;


    struct vdescriptor
    {
      using type = siconos::any::vdescriptor;
    };
    collection<vdescriptor> _vdescriptors;

    std::array<typename env::template collection<item_types>, 1> _items;

    /* methods & operators */
    template<typename Fun>
    decltype (auto) operator () (Fun&& fun)
    {
      return proj(*this)(fun);
    }

  };

  template<typename Env, typename OSI, typename ...Items>
  static constexpr auto make_data()
  {
    return data<Env,OSI, Items...>{};
  }

  static constexpr void for_each(auto&& fun, auto&& tpl)
  {
    std::apply([&fun](auto&&... args) { ((fun(args)), ...);}, tpl);
  };

  template<typename T>
  static constexpr void for_each_attribute(auto&& fun, auto& data)
  {
    for_each(
      [&fun,&data]<typename ...Attrs>(Attrs&...)
      {
        (fun(siconos::get_array<Attrs>(data)), ...);
      },
      typename T::attributes{});
  };

  template<typename T>
  auto add_vertex_item(const auto step, auto& data)
  {
    auto vd = data._graph.add_vertex(data._counter++);

    data._vdescriptors[step%std::size(data._vdescriptors)].push_back(vd);
    data._items[step%std::size(data._items)].push_back(T{});

    data._graph.index(vd) = std::size(data._vdescriptors[step%std::size(data._vdescriptors)])-1;

    for_each_attribute<T>(
      [&step](auto& a)
      {
        a[step%std::size(a)].push_back(typename std::decay_t<decltype(a)>::value_type::value_type{});
      },
      data);

    return vd;
  }

  void remove_vertex_item(const auto vd,
                          const auto step,
                          auto& data)
  {
    auto i = data._graph.index(vd);

    data._graph.remove_vertex(vd);

    typename std::decay_t<decltype(data)>::item_types item =
      data._items[step%std::size(data._items)][i];

    move_back(i, data._vdescriptors[step%std::size(data._vdescriptors)]);
    move_back(i, data._items[step%std::size(data._items)]);
    data._graph.index(data._vdescriptors[step%std::size(data._vdescriptors)][i]) = i;

    std::visit([&i,&step,&data]<typename Item>(Item)
    {
      for_each_attribute<Item>(
        [&i,&step](auto &a)
        {
          move_back(i, a[step%std::size(a)]);
        },
        data);
    }, item);
  }

  template<typename T>
  void add_item(const auto step, auto& data)
  {
    for_each_attribute<T>(
      [&step](auto& a)
      {
        a[step%std::size(a)].push_back(typename std::decay_t<decltype(a)>::value_type::value_type{});
      },
      data);
  }
};


#endif
