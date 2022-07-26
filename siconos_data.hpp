#ifndef SICONOS_DATA_HPP
#define SICONOS_DATA_HPP

#include <type_traits>
#include <variant>
#include <concepts>

#include "siconos_util.hpp"
#include "siconos_traits.hpp"
namespace siconos
{
  template<typename T, std::size_t N>
  using memory_t = std::array<T, N>;

  template<typename Tag>
  static constexpr auto memory_size = [](concepts::tuple_like auto&& keeps)
    constexpr -> decltype(auto)
  {
    // let rec loop(...) = ...
    constexpr auto rec_loop = []<typename Tpl>(auto&& loop, Tpl) constexpr
    {
      constexpr auto keep = car(Tpl{});
      using keep_t = std::decay_t<decltype(keep)>;

      if constexpr (std::is_same_v<Tag, typename keep_t::tag>)
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

    // at least one element remaining
    static_assert(std::tuple_size_v<std::decay_t<decltype(keeps)>> > 0);

    return rec_loop(rec_loop, keeps);
  };

  static constexpr auto memory =
    []<typename T>(typename T::size_type step, T& mem)
    constexpr -> typename T::value_type&
  {
    return mem[step%std::size(mem)];
  };

  template<concepts::tag T>
  static constexpr auto get_memory = [](auto& data) constexpr
    -> decltype(auto)
  {
    using data_t = std::decay_t<decltype(data)>;
    constexpr typename data_t::data_types t = T{};
    auto& mem = std::get<t.index()>(data._collections);
    return mem;
  };

  static constexpr auto get_step = [](auto& data) constexpr -> decltype(auto)
  {
    using data_t = std::decay_t<decltype(data)>;
    using time_discretization = typename data_t::time_discretization;
    // one time_discretization / data
    return get_memory<typename time_discretization::step>(data)[0][0];
  };

  template<typename T>
  static constexpr auto get_step_data = [](const auto vd, const auto step, auto& data) -> decltype(auto)
  {
    auto indx = data._graph.index(vd);
    auto& mem = get_memory<T>(data);
    return memory(step, mem)[indx];
  };

  template<typename T>
  static constexpr auto get = [](const auto vd, auto& data) -> decltype(auto)
  {
    auto step = get_step(data);
    auto indx = data._graph.index(vd);
    auto& mem = get_memory<T>(data);
    return memory(step, mem)[indx];
  };

  static constexpr auto move_back = [](const auto i, auto& a)
  {
    a[i] = std::move(a.back());
    a.pop_back();
  };

  template<typename Env, typename Sim, typename ...Attributes>
  struct data
  {
    using env = Env;
    using osi = typename Sim::one_step_integrator;
    using time_discretization = typename Sim::time_discretization;

    using indice = typename env::indice;
    using scalar = typename env::scalar;
    using graph_t = typename env::graph;

    // wrap collection type in a memory container according to
    // osi::keeps specification
    template<typename C>
    using collection =
      memory_t<typename env::template collection<
                 typename traits::config<
                   env, typename C::structure::type>::type>,
               memory_size<typename C::tag>(typename osi::definition::keeps{})>;

    // global variant for all attributes
    using data_types = std::variant<typename Attributes::tag ...>;

    // global storage type for all defined collections
    using collections_t = std::tuple<collection<Attributes>...>;

    // global storage
    collections_t _collections;

    // administrative counter
    indice _counter = 0;

    // graph
    graph_t _graph;

    using vertex_items =
      decltype(
        std::apply(
          []<typename ...Vertex_items>(Vertex_items...)
          {
            return std::variant<Vertex_items ...>{};
          }, all_vertex_items(Sim{})));

    // vertices descriptors
    struct vdescriptor : some::tag {};
    struct vdescriptor_attr :
      attribute<
      tag<vdescriptor>,
      symbol<"vd">,
      description<"vertex descriptor">,
      structure<some::vdescriptor<vertex_items>>> {};

    collection<vdescriptor_attr> _vdescriptors;

    ///////////////////////
    // methods & operators
    template<typename Fun>
    decltype (auto) operator () (Fun&& fun)
    {
      return proj(*this)(fun);
    }

  };


  template<typename Env, typename Sim>
  static auto make_data = []() constexpr
  {
    using attributes = decltype(all_terminal_attributes(Sim{}));

    return std::apply(
      []<typename ...Attributes>(Attributes...)
      {
        return data<Env, Sim, Attributes...>{};
      }, attributes{});
  };

  template<concepts::item T>
  static auto for_each_attribute = [](auto&& fun, auto& data)
    constexpr
  {
    for_each(
      [&fun,&data]<concepts::attribute ...Attrs>(Attrs&...)
      {
        (fun(siconos::get_memory<typename Attrs::tag>(data)), ...);
      },
      all_terminal_attributes(T{}));
  };

  template<typename T>
  static auto add_attributes = [](auto& data) constexpr -> decltype(auto)
  {
    auto step = get_step(data);

    for_each_attribute<T>(
      [&step](auto& a)
      {
        memory(step, a).push_back(typename std::decay_t<decltype(a)>::value_type::value_type{});
      },
      data);

  };

  template<typename T>
  static auto add_vertex_item = [](auto& data) constexpr -> decltype(auto)
  {
    static_assert(concepts::vertex_item_t<T, std::decay_t<decltype(data)>>);

    auto step = get_step(data);

    auto vd = data._graph.add_vertex(data._counter++);

    memory(step, data._vdescriptors).push_back({vd, T{}});
//    memory(step, data._items).push_back(T{});

    data._graph.index(vd) = std::size(memory(step, data._vdescriptors))-1;

    add_attributes<T>(data);
    return vd;
  };

  template<typename T>
  static auto add = [](auto& data) constexpr -> decltype(auto)
  {
    if constexpr (concepts::vertex_item_t<T, std::decay_t<decltype(data)>>)
    {
      return add_vertex_item<T>(data);
    }
    else
    {
      add_attributes<T>(data);
    }
  };

  static auto remove_vertex_item = [](const auto vd,
                                      const auto step,
                                      auto& data) constexpr
  {

    auto i = data._graph.index(vd);

    data._graph.remove_vertex(vd);

    move_back(i, memory(step, data._vdescriptors));
    data._graph.index(std::get<0>((memory(step, data._vdescriptors)[i]).value)) = i;

    std::visit([&i,&step,&data]<typename Item>(Item)
    {
      for_each_attribute<Item>(
        [&i,&step](auto &a)
        {
          move_back(i, memory(step, a));
        },
        data);
    }, (std::get<1>((memory(step, data._vdescriptors)[i]).value)));
  };

  template<typename T>
  static auto build_a_value = [](auto& data) constexpr -> decltype(auto)
  {
    using data_t = std::decay_t<decltype(data)>;
    constexpr typename data_t::data_types t = T{};
    return t;
  };

};

#endif
