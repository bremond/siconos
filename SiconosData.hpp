#ifndef SICONOS_DATA_HPP
#define SICONOS_DATA_HPP

#include <variant>
#include <tuple>
#include <boost/pfr.hpp>
#include <boost/hana.hpp>
#include <boost/hana/ext/std/tuple.hpp>

namespace siconos
{

  template<typename OSI, typename ...Ds>
  struct data
  {
    using env = typename OSI::Env;
    using indice = env::indice;
    using scalar = env::scalar;

    template<typename T>
    using collection = std::array<typename env::collection<T>,
                                  typename OSI::memory_size>;

    using types = std::variant<Ds...>;
    std::tuple<collection<typename Ds::type>...> collections;

    env::indice counter = 0;
    env::graph graph;

    collection<typename env::vdescriptor> vdescriptors;

  };

  template<typename OSI>
  static constexpr auto make_data()
  {
    using env = typename OSI::env;
    using model = OSI::model;
    return std::apply([]<typename ...Ds>(Ds&&...) { return data<OSI,Ds...>{}; },
                      boost::pfr::structure_to_tuple(model{}));

  }

  template<typename T, typename D>
  typename D::collection<typename T::type>::value_type&
  get(const T, D& data)
  {
    constexpr typename D::types t = T{};
    return std::get<t.index()>(data.collections).back();
  };



  void for_each(auto&& fun, auto& collections)
  {
    std::apply([&fun](auto&&... args) { ((fun(args)), ...);}, collections);
  };

  auto add_item(auto& data)
  {
    auto vd = data.graph.add_vertex(data.counter++);
    data.vdescriptors.push_back(vd);

    for_each(
      [](auto& a)
      { a.push_back(typename std::decay_t<decltype(a)>::value_type{}); },
      data.collections);

    return std::size(data.vdescriptors);
  }

  void move_back(const auto i, auto& a)
  {
    a[i] = std::move(a.back());
    a.pop_back();
  }

  void remove_item(const auto i, auto& data)
  {
    auto& vd = data.vdescriptors[i];

    data.graph.remove_vertex(vd);

    move_back(i, data.vdescriptors);
    for_each([&i](auto &a) { move_back(i, a);}, data.collections);
  }
};


#endif
