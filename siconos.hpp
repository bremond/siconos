#ifndef SICONOS_HPP
#define SICONOS_HPP

#include "siconos_storage.hpp"
#include "siconos_pattern.hpp"
#include "siconos_linear_algebra.hpp"
#include <range/v3/view/zip.hpp>
namespace siconos
{
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

      // should not be an attribute
      struct fext : some::vector<dof>, // some::function<...>
                    access<fext> {};

      using attributes = gather<mass_matrix, q, velocity, fext>;

    };

    struct relation : item<>
    {
      struct h_matrix : some::matrix<1, dof>, access<h_matrix> {};

      using attributes = gather<h_matrix>;

      static constexpr auto compute_input = []()
      {
      };

      static constexpr auto compute_output = []<match::item Interaction,
                                                std::size_t deriv_num>
        (double time, Interaction, auto& data)
      {
        auto& y = get<typename Interaction::y>(data);
        auto& H = get<h_matrix>(data);
      };

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

  template<typename Formulation, typename Interaction>
  struct topology : item<>
  {
    using dynamical_system = typename Formulation::dynamical_system;
    using interaction = Interaction;

    struct dsg0 : some::graph<some::item_ref<dynamical_system>,
                              some::item_ref<interaction>>{};

    struct ig0 : some::graph<some::item_ref<interaction>,
                             some::item_ref<dynamical_system>>{};

    struct ig1 : some::graph<some::item_ref<interaction>,
                             some::item_ref<dynamical_system>>{};

    using attributes = gather<dsg0, ig0, ig1>;
  };

  template<typename Nslaw, typename Formulation, std::size_t N>
  struct interaction : item<>
  {
    struct nonsmooth_law : some::item_ref<Nslaw>, access<nonsmooth_law> {};
    struct relation      : some::item_ref<typename Formulation::relation>,
                           access<relation> {};

    struct y : some::vector<N*Formulation::dof>{};
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

      };

    };

    struct moreau_jean : item<>
    {
      struct theta : some::scalar {};

      using attributes = gather<theta>;

      using kinds = gather<
        keep<typename system::q, 2>,
        keep<typename system::velocity, 2>>;

      static constexpr auto compute_iteration_matrix = [](auto step, auto& data)
        constexpr -> decltype(auto)
      {
        auto& mass_matrices = get_memory<mass_matrix>(data);
        auto& external_forces = get_memory<fext>(data);

        auto& Ms =      memory(step, mass_matrices);
        auto& fs =      memory(step, external_forces);

        if constexpr(is_kind_of<mass_matrix, some::time_invariant>(data))
        {
          if constexpr(is_kind_of<fext, some::time_invariant>(data))
          {
            if constexpr(is_kind_of<mass_matrix, some::diagonal>(data))
            {
              for (auto [M, f] : ranges::views::zip(Ms, fs))
              {
                linear_algebra::solve_in_place(M, f);
              }
            }
          }
        }
      };

      static void compute_free_state(auto step, auto h, auto& data)
      {
         auto& velocities = get_memory<velocity>(data);
         auto fexts  = get_memory<fext>(data);
         auto theta = get<moreau_jean::theta>(data);

         auto& vs =      memory(step, velocities);
         auto& vs_next = memory(step+1, velocities);
         auto& minv_fs      = memory(step, fexts);
         auto& minv_fs_next      = memory(step+1, fexts);

         // f <- M  f
         for (auto [v, v_next, minv_f, minv_f_next] :
                ranges::views::zip(vs, vs_next, minv_fs, minv_fs_next))
         {
           v_next = v + h*theta*minv_f + h*(1 - theta)*minv_f_next;
         }
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

  template<typename TD, typename OSI, typename OSNSPB, typename TOPO>
  struct time_stepping : item<>
  {
    using time_discretization = TD;
    using one_step_integrator = OSI;
    using one_step_nonsmooth_problem = OSNSPB;
    using topology = TOPO;

    using items = gather<time_discretization,
                         one_step_integrator,
                         one_step_nonsmooth_problem,
                         topology>;

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
