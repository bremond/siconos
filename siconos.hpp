#ifndef SICONOS_HPP
#define SICONOS_HPP

#include <vector>
#include <array>
#include "SiconosGraph.hpp"

#include "siconos_storage.hpp"
#include "siconos_pattern.hpp"

namespace siconos
{
  template<typename Param>
  struct lagrangian
  {
    static constexpr auto dof = Param::dof;

    struct dynamical_system : vertex_item<
      description
      <"the dynamical system">>
    {
      struct mass_matrix :
        some::matrix<dof, dof> {};

      struct q : some::vector<dof> {};

      struct velocity : some::vector<dof> {};

      struct fext : some::vector<dof> {};

      using attributes = gather<mass_matrix, q, velocity, fext>;

    };

    struct relation : item<>
    {
      struct h_matrix : some::matrix<1, dof> {};

      using attributes = gather<h_matrix>;
    };
  };

  struct nonsmooth_law
  {
    struct newton_impact_friction : item<>
    {
      struct e : some::scalar {};
      struct mu : some::scalar {};

      using attributes = gather<e, mu>;
    };

    struct newton_impact : item<>
    {
      struct e : some::scalar {};

      using attributes = gather<e>;

    };
  };

  template<typename Nslaw, typename Relation>
  struct interaction : edge_item<>
  {
    using nonsmooth_law = Nslaw;
    using relation = Relation;

    using uses = gather<nonsmooth_law,
                        relation>;
  };

  struct lcp
  {};

  template<typename Type>
  struct one_step_nonsmooth_problem : item<>
  {
    using problem_type = Type;

    struct level : some::indice {};

    using attributes = gather<level>;
  };

  template<typename Form>
  struct one_step_integrator
  {
    using formulation = Form;
    using system = typename formulation::dynamical_system;

    struct moreau_jean : item<>
    {
      struct theta : some::scalar {};

      using attributes = gather<theta>;

      using keeps = gather<
        keep<typename system::q, 2>,
        keep<typename system::velocity, 2>>;
    };
  };

  template<typename Param>
  struct time_discretization : item<>
  {
    struct step : some::indice {};
    using attributes = tuple<step>;
  };

  template<typename TD, typename OSI, typename OSNSPB, typename ...GraphItems>
  struct time_stepping : item<>
  {
    using time_discretization = TD;
    using one_step_integrator = OSI;
    using one_step_nonsmooth_problem = OSNSPB;
    using graph_items = gather<GraphItems...>;

    using uses = gather<time_discretization,
                        one_step_integrator,
                        one_step_nonsmooth_problem,
                        GraphItems...>;
  };

}

#endif
