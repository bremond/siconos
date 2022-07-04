#ifndef SICONOS_HPP
#define SICONOS_HPP

#include <vector>
#include <array>
#include "SiconosGraph.hpp"

#include "siconos_data.hpp"

namespace siconos
{
  template<typename X>
  struct item
  {
    using type = X;

    static constexpr auto add(const auto step, auto& data)
    {
      return siconos::add_vertex_item(step, data);
    };
  };


  template<typename T>
  struct attribute
  {
    using type = T;

    static constexpr auto get = [](const auto vd, const auto step, auto& data) -> typename T::type&
    {
//      return siconos::get<T>(vd, step, data);
    };
  };

}

struct env
{
  using scalar = double;
  using indice = std::size_t;

  using graph = SiconosGraph < indice, indice,
                               boost::no_property,
                               boost::no_property,
                               boost::no_property >;

  using vdescriptor = graph::VDescriptor;

  template<typename T>
  using collection = std::vector<T>;

  template<indice N>
  using vector = std::array<scalar, N>;

  template<indice N, indice M>
  using matrix = std::array<scalar, N*M>;

};

struct param
{
  static constexpr auto dof = 3;
};

namespace siconos
{
  template<typename ...Args>
  using tuple = std::tuple<Args...>;

  template<typename Param>
  struct lagrangian
  {
    static constexpr auto dof = Param::dof;

    struct dynamical_system
    {
      struct mass_matrix
      {
        using type = any::matrix<dof, dof>;
      };

      struct q
      {
        using type = any::vector<dof>;
      };

      struct velocity
      {
        using type = any::vector<dof>;
      };

      struct fext
      {
        using type = any::vector<dof>;
      };

      using attributes =
        tuple<mass_matrix, q, velocity, fext>;

    };

    struct relation
    {
    };

    using vertex_items = tuple<dynamical_system>;
    using edge_items = tuple<relation>;
  };

  template<typename NSLAW, typename RELATION>
  struct interaction
  {
    using nonsmooth_law = NSLAW;
    using relation = RELATION;
  };

  struct nonsmooth_law
  {

    struct newton_impact_friction
    {
      struct e
      {
        using type = any::scalar;
      };

      struct mu
      {
        using type = any::scalar;
      };

      using attributes = tuple<e, mu>;
    };

    struct newton_impact
    {
      struct e
      {
        using type = any::scalar;
      };

      using attributes = tuple<e>;
    };

    using items = tuple<newton_impact_friction, newton_impact>;
  };

  struct lcp
  {};

  template<typename Type>
  struct one_step_nonsmooth_problem
  {
    using type = Type;
  };

  template<typename Form>
  struct one_step_integrator
  {
    struct moreau_jean
    {
      using formulation = Form;
      using system = typename formulation::dynamical_system;

      template<std::size_t N, typename ...As>
      struct keeper
      {
        static constexpr std::size_t value = N;
        using type = tuple<As...>;
      };

      using keep = keeper<2, typename system::q, typename system::velocity>;

      struct theta
      {
        using type = any::scalar;
      };

      using attributes = tuple<theta>;
    };

    using items = tuple<moreau_jean>;
  };

  template<typename Param>
  struct time_discretization
  {
//    struct event_manager
//    {
      struct current_time_step
      {
        using type = any::indice;
      };
//    };

    using attributes = tuple<current_time_step>;
  };

  template<typename TD, typename OSI, typename OSNSPB>
  struct time_stepping
  {
    using time_discretization = TD;
    using one_step_integrator = OSI;
    using one_step_nonsmooth_problem = OSNSPB;

    using attributes = tuple<typename time_discretization::current_time_step>;
  };
}

#endif
