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
  template<typename OSI, typename ...Ds>
  struct data
  {
    using osi = OSI;
    using env = osi::env;
    using indice = env::indice;
    using scalar = env::scalar;
    using graph_t = env::graph;

    template<typename T>
    using collection =
      std::conditional<contains<T>(typename osi::keep::type{}),
                       std::array<typename env::collection<T>, osi::keep::value>,
                       std::array<typename env::collection<T>, 1>>::type;

    using types = std::variant<Ds...>;
    std::tuple<collection<typename Ds::type>...> collections;

    indice counter = 0;
    graph_t graph;

    collection<typename env::vdescriptor> vdescriptors;

  };

  template<typename OSI>
  static constexpr auto make_data()
  {
    using formulation = OSI::formulation;
    using env = OSI::env;

    return std::apply([]<typename ...Ds>(Ds&&...) { return data<OSI, Ds...>{}; },
                      typename formulation::dynamical_system<env>::attributes{});

  }

  template<typename T, typename D>
  auto&
  get(const typename D::env::indice step, D& data)
  {
    constexpr typename D::types t = T{};
    auto& array = std::get<t.index()>(data.collections);
    return array[step % std::size(array)].back();
  };

  void for_each(auto&& fun, auto& collections)
  {
    std::apply([&fun](auto&&... args) { ((fun(args)), ...);}, collections);
  };

  auto add_item(const auto step, auto& data)
  {
    auto vd = data.graph.add_vertex(data.counter++);
    data.vdescriptors[step%std::size(data.vdescriptors)].push_back(vd);

    for_each(
      [&step](auto& a)
      { a[step%std::size(a)].push_back(typename std::decay_t<decltype(a)>::value_type::value_type{}); },
      data.collections);

    return std::size(data.vdescriptors);
  }

  void move_back(const auto i, auto& a)
  {
    a[i] = std::move(a.back());
    a.pop_back();
  }

  void remove_item(const auto i, const auto step, auto& data)
  {
    auto& vd = data.vdescriptors[step%std::size(data.vdescriptors)][i];

    data.graph.remove_vertex(vd);

    move_back(i, data.vdescriptors[step%std::size(data.vdescriptors)]);
    for_each([&i,&step](auto &a) { move_back(i, a[step%std::size(a)]);}, data.collections);
  }
};


#endif
