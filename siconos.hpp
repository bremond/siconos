#ifndef SICONOS_HPP
#define SICONOS_HPP

#include <vector>
#include <array>
#include "SiconosGraph.hpp"

#include "siconos_data.hpp"
#include "siconos_pattern.hpp"

namespace siconos
{
  template<typename Param>
  struct lagrangian
  {
    static constexpr auto dof = Param::dof;

    struct dynamical_system
    {
      struct mass_matrix : some::tag {};
      struct q : some::tag {};
      struct velocity : some::tag {};
      struct fext : some::tag {};

      using definition = vertex_item<
        description
        <"the dynamical system">,

        attribute<
          tag<mass_matrix>,
          symbol<"M">,
          description<"the mass matrix">,
          structure<some::matrix<dof, dof>>>,

        attribute<
          tag<q>,
          symbol<"q">,
          description<"state">,
          structure<some::vector<dof>>>,

        attribute<
          tag<velocity>,
          symbol<"v">,
          description<"the velocity">,
          structure<some::vector<dof>>>,

        attribute<
          tag<fext>,
          symbol<"fext">,
          description<"external force">,
          structure<some::vector<dof>>>>;
    };

    struct relation
    {
      struct h_matrix : some::tag {};

      using definition = item<

        description
        <"the relation">,

        attribute<
          tag<h_matrix>,
          symbol<"H">,
          description<"the H matrix">,
          structure<some::matrix<1, dof>>>>;
    };
  };

  struct nonsmooth_law
  {
    struct newton_impact_friction
    {
      struct e : some::tag {};
      struct mu : some::tag {};

      using definition = item<

        description
        <"the Newton impact friction law">,

        attribute<
          tag<e>,
          symbol<"e">,
          description<"restitution coefficient">,
          structure<some::scalar>>,

        attribute<
          tag<mu>,
          symbol<"mu">,
          description<"Coulomb friction coefficient">,
          structure<some::scalar>>>;
    };

    struct newton_impact
    {
      struct e : some::tag {};

      using definition = item<

        description
        <"The Newton impact law">,

        attribute<
          tag<e>,
          symbol<"e">,
          description<"restitution coefficient">,
          structure<some::scalar>>>;
    };
  };

  template<typename Nslaw, typename Relation>
  struct interaction
  {
    using nonsmooth_law = Nslaw;
    using relation = Relation;

    using definition = item<
      description
      <"The interaction">,

      use<nonsmooth_law>,
      use<relation>>;
  };

  struct lcp
  {};

  template<typename Type>
  struct one_step_nonsmooth_problem
  {
    using problem_type = Type;

    struct level : some::tag {};

    using definition = item<

      description
      <"The one step nonsmooth problem">,
      attribute<
        tag<level>,
        symbol<"level">,
        description<"Index set level">,
        structure<some::indice>>>;
  };

  template<typename Form>
  struct one_step_integrator
  {
    using formulation = Form;
    using system = typename formulation::dynamical_system;

    struct moreau_jean
    {
      struct theta : some::tag {};

      using definition = item<

        description
        <"The Moreau-Jean integrator">,

        use<system>,
        attribute<
          tag<theta>,
          symbol<"theta">,
          description<"theta method parameter">,
          structure<some::scalar>>,

        keep<typename system::q, 2>,
        keep<typename system::velocity, 2>>;
    };
  };

  template<typename Param>
  struct time_discretization
  {
    struct step : some::tag {};
    using definition = item<

      description
      <"the time discretization">,

        attribute<
          tag<step>,
          symbol<"step">,
          description<"">,
          structure<some::indice>>>;

  };

  template<typename TD, typename OSI, typename OSNSPB, typename ...Inters>
  struct time_stepping
  {
    using time_discretization = TD;
    using one_step_integrator = OSI;
    using one_step_nonsmooth_problem = OSNSPB;
    using interactions = std::tuple<Inters...>;

    using definition = item<
      description
      <"The time stepping">,
      use<time_discretization>,
      use<one_step_integrator>,
      use<one_step_nonsmooth_problem>,
      use<Inters...>>;
  };
}

#endif
