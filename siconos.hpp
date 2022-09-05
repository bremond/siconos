#ifndef SICONOS_HPP
#define SICONOS_HPP

#include "siconos_storage.hpp"
#include "siconos_pattern.hpp"

namespace siconos
{
  template<typename ...Params>
  struct lagrangian : frame<Params...>
  {
    static constexpr auto dof = frame<Params...>::dof;

    struct dynamical_system : vertex_item<
      description
      <"A lagrangian dynamical system [...]">>
    {
      struct mass_matrix :
        some::matrix<dof, dof> {};

      struct q : some::vector<dof>,
                 access<q> {};

      struct velocity : some::vector<dof>,
                        access<velocity> {};

      struct fext : some::vector<dof>,
                    access<fext> {};

      using attributes = gather<mass_matrix, q, velocity, fext>;

    };

    struct relation : item<>
    {
      struct h_matrix : some::matrix<1, dof>, access<h_matrix> {};

      using attributes = gather<h_matrix>;
    };
  };

  struct nonsmooth_law
  {
    struct newton_impact_friction : item<>
    {
      struct e : access<e>, some::scalar {};
      struct mu : some::scalar, access<mu> {};

      using attributes = gather<e, mu>;
    };

    struct newton_impact : item<>
    {
      struct e : some::scalar, access<e> {};

      using attributes = gather<e>;

    };
  };

  template<typename Nslaw, typename Relation>
  struct interaction : edge_item<>
  {
    struct nonsmooth_law : some::item_ref<Nslaw>, access<nonsmooth_law> {};
    struct relation : some::item_ref<Relation>, access<relation> {};

    using attributes = gather<nonsmooth_law,
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

  template<typename ...Params>
  struct time_discretization : item<>
  {
    struct t0 : some::scalar {};
    struct h : some::scalar {};
    struct step : some::indice {};
    using attributes = gather<h, t0, step>;
  };

  template<typename TD, typename OSI, typename OSNSPB, typename ...GraphItems>
  struct time_stepping : item<>
  {
    using time_discretization = TD;
    using one_step_integrator = OSI;
    using one_step_nonsmooth_problem = OSNSPB;
    using graph_items = gather<GraphItems...>;

    using items = gather<time_discretization,
                         one_step_integrator,
                         one_step_nonsmooth_problem,
                         GraphItems...>;
  };

}

#endif
