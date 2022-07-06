#ifndef SICONOS_DATA_HPP
#define SICONOS_DATA_HPP

#include <type_traits>
#include <variant>
#include <tuple>

namespace siconos
{
  struct any
  {
    struct scalar{};
    struct indice{};
    struct vdescriptor{};

    template<std::size_t N, std::size_t M>
    struct matrix{};
    template<std::size_t N>
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

  template<typename T, std::size_t N>
  using memory_t = std::array<T, N>;

  template<typename T>
  typename T::value_type& memory(auto step, T& mem)
  {
    return mem[step%std::size(mem)];
  };

  template<typename T>
  constexpr auto get_memory = [](auto& data) -> decltype(auto)
  {
    using data_t = std::decay_t<decltype(data)>;
    constexpr typename data_t::data_types t = T{};
    auto& mem = std::get<t.index()>(data._collections);
    return mem;
  };

  template<typename T>
  constexpr auto get = [](const auto vd, const auto step, auto& data) -> decltype(auto)
  {
    auto indx = data._graph.index(vd);
    auto& mem = get_memory<T>(data);
    return memory(step, mem)[indx];
  };

  template<typename T>
  constexpr auto xget = [](const auto vd, auto& data) -> decltype(auto)
  {
    using data_t = std::decay_t<decltype(data)>;
    using time_discretization = typename data_t::time_discretization;
    auto step = get_memory<typename time_discretization::current_time_step>(data)[0][0];
    auto indx = data._graph.index(vd);
    auto& mem = get_memory<T>(data);
    return memory(step, mem)[indx];
  };

  void move_back(const auto i, auto& a)
  {
    a[i] = std::move(a.back());
    a.pop_back();
  }

  template<typename Env, typename Sim, typename ...Items>
  struct data
  {
    using env = Env;
    using osi = typename Sim::one_step_integrator;
    using time_discretization = typename Sim::time_discretization;

    using indice = typename env::indice;
    using scalar = typename env::scalar;
    using graph_t = typename env::graph;

    template<typename C>
    using collection =
      typename std::conditional<contains<C>(typename osi::keep::type{}),
                                memory_t<typename env::template collection<typename traits::config<env, typename C::type>::type>, osi::keep::value>,
                                memory_t<typename env::template collection<typename traits::config<env, typename C::type>::type>, 1>>::type;

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

    memory_t<typename env::template collection<item_types>, 1> _items;

    /* methods & operators */
    template<typename Fun>
    decltype (auto) operator () (Fun&& fun)
    {
      return proj(*this)(fun);
    }

  };

  template<typename Env, typename Sim, typename Inter>
  static constexpr auto make_data()
  {
    using system = typename Sim::one_step_integrator::system;
    using nslaw = typename Inter::nonsmooth_law;
    using time_discretization = typename Sim::time_discretization;
    return data<Env, Sim, system, nslaw, time_discretization>{};
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
        (fun(siconos::get_memory<Attrs>(data)), ...);
      },
      typename T::attributes{});
  };

  template<typename T>
  auto add_vertex_item(const auto step, auto& data)
  {
    auto vd = data._graph.add_vertex(data._counter++);

    memory(step, data._vdescriptors).push_back(vd);
    memory(step, data._items).push_back(T{});

    data._graph.index(vd) = std::size(memory(step, data._vdescriptors))-1;

    for_each_attribute<T>(
      [&step](auto& a)
      {
        memory(step, a).push_back(typename std::decay_t<decltype(a)>::value_type::value_type{});
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
      memory(step, data._items)[i];

    move_back(i, memory(step, data._vdescriptors));
    move_back(i, memory(step, data._items));
    data._graph.index(memory(step, data._vdescriptors)[i]) = i;

    std::visit([&i,&step,&data]<typename Item>(Item)
    {
      for_each_attribute<Item>(
        [&i,&step](auto &a)
        {
          move_back(i, memory(step, a));
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
        memory(step, a).push_back(typename std::decay_t<decltype(a)>::value_type::value_type{});
      },
      data);
  }
};

#endif
