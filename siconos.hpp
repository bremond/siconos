#ifndef SICONOS_HPP
#define SICONOS_HPP

#include "siconos_storage.hpp"
#include "siconos_pattern.hpp"
#include "siconos_linear_algebra.hpp"
#include <range/v3/view/zip.hpp>
namespace siconos
{
  struct graph_item {};

  template<typename ...Params>
  struct lagrangian : frame<Params ...>
  {
    static constexpr auto dof = frame<Params...>::dof;

    struct dynamical_system : item<
      description
      <"A lagrangian dynamical system [...]">>
    {
      struct mass_matrix :
        some::matrix<dof, dof>,
        access<mass_matrix> {};

      struct q : some::vector<dof>,
                 access<q> {};

      struct velocity : some::vector<dof>,
                        access<velocity> {};

      struct fext : some::vector<dof>, // some::function<...>
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

  struct topology : item<>
  {
    struct dsg0 : some::graph {};
    struct ig0 : some::graph {};
    struct ig1 : some::graph {};

    using attributes = gather<dsg0, ig0, ig1>;
  };

  template<typename Nslaw, typename Relation>
  struct interaction : item<>
  {
    struct nonsmooth_law : some::item_ref<Nslaw>, access<nonsmooth_law> {};
    struct relation      : some::item_ref<Relation>, access<relation> {};

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
    using q = typename system::q;
    using velocity = typename system::velocity;
    using mass_matrix = typename system::mass_matrix;
    using fext = typename system::fext;

    struct euler : item<>
    {
      using kinds = gather<
        keep<q, 2>,
        keep<velocity, 2>>;

      static void compute_free_state(auto step, auto h, auto& data)
      {
        auto& velocities = get_memory<velocity>(data);
        auto& mass_matrices = get_memory<mass_matrix>(data);
        auto& external_forces = get_memory<fext>(data);

        auto& Ms =      memory(step, mass_matrices);
        auto& vs =      memory(step, velocities);
        auto& vs_next = memory(step+1, velocities);
        auto& fs =      memory(step, external_forces);

        for (auto [M, v, v_next, f] : ranges::views::zip(Ms, vs, vs_next, fs))
        {
          v_next = v + h * linear_algebra::solve(M, f);
        }

      };

    };

    struct moreau_jean : item<>
    {
      struct theta : some::scalar {};

      using attributes = gather<theta>;

      using kinds = gather<
        keep<typename system::q, 2>,
        keep<typename system::velocity, 2>>;

      static void compute_free_state(auto step, auto h, auto& data)
      {
      };
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

  template<typename TD, typename OSI, typename OSNSPB, typename TOPO, typename ...GraphItems>
  struct time_stepping : item<>
  {
    using time_discretization = TD;
    using one_step_integrator = OSI;
    using one_step_nonsmooth_problem = OSNSPB;
    using topology = TOPO;
    using graph_items = gather<GraphItems...>;

    using items = gather<time_discretization,
                         one_step_integrator,
                         one_step_nonsmooth_problem,
                         topology,
                         GraphItems...>;

    static decltype(auto) current_step(auto& data)
    {
      return get<typename time_discretization::step>(data)[0];
    }

    static decltype(auto) time_step(auto& data)
    {
      return get<typename time_discretization::h>(data)[0];
    }

    static void compute_one_step(auto& data)
    {
      one_step_integrator::compute_free_state(current_step(data),
                                              time_step(data), data);
      current_step(data) += 1;

    }
  };

}

#endif
