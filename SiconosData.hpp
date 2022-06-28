#ifndef SICONOS_DATA_HPP
#define SICONOS_DATA_HPP

#include <type_traits>
#include <variant>
#include <tuple>
#include <boost/pfr.hpp>
#include <boost/hana.hpp>
#include <boost/hana/ext/std/tuple.hpp>

template<typename U, typename... T>
constexpr bool contains(std::tuple<T...>)
{
  return (std::is_same_v<U, T> || ...);
}

namespace siconos
{
  struct types
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
    struct config<E, types::scalar>
    {
      using type = typename E::scalar;
    };

    template<typename E>
    struct config<E, types::indice>
    {
      using type = typename E::indice;
    };

    template<typename E>
    struct config<E, types::vdescriptor>
    {
      using type = typename E::vdescriptor;
    };

    template<std::size_t N, typename E>
    struct config<E, types::vector<N>>
    {
      using type = typename E::template vector<N>;
    };

    template<std::size_t N, std::size_t M, typename E>
    struct config<E, types::matrix<N, M>>
    {
      using type = typename E::template matrix<N, M>;
    };
  };

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

  template<typename T>
  auto get = [](const auto vd, const auto step, auto& data) -> decltype(auto)
  {
    using data_t = std::decay_t<decltype(data)>;
    constexpr typename data_t::data_types t = T{};
    auto indx = data._graph.index(vd);
    auto& array = std::get<t.index()>(data._collections);
    return array[step % std::size(array)][indx];
  };

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
    collections_t _collections;

    indice _counter = 0;
    graph_t _graph;

    struct vdescriptor
    {
      using type = siconos::types::vdescriptor;
    };
    collection<vdescriptor> vdescriptors;


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


  void for_each(auto&& fun, auto& collections)
  {
    std::apply([&fun](auto&&... args) { ((fun(args)), ...);}, collections);
  };

  auto add_item(const auto step, auto& data)
  {
    auto vd = data._graph.add_vertex(data._counter++);
    data.vdescriptors[step%std::size(data.vdescriptors)].push_back(vd);
    data._graph.index(vd) = std::size(data.vdescriptors[step%std::size(data.vdescriptors)])-1;

    for_each(
      [&step](auto& a)
      { a[step%std::size(a)].push_back(typename std::decay_t<decltype(a)>::value_type::value_type{}); },
      data._collections);

    return vd;
  }

  void move_back(const auto i, auto& a)
  {
    a[i] = std::move(a.back());
    a.pop_back();
  }

  void remove_item(const auto i, const auto step, auto& data)
  {
    auto& vd = data.vdescriptors[step%std::size(data.vdescriptors)][i];

    data._graph.remove_vertex(vd);

    move_back(i, data.vdescriptors[step%std::size(data.vdescriptors)]);
    data._graph.index(data.vdescriptors[step%std::size(data.vdescriptors)][i]) = i;

    for_each([&i,&step](auto &a) { move_back(i, a[step%std::size(a)]);}, data._collections);
  }
};


#endif
