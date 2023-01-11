#ifndef SICONOS_HPP
#define SICONOS_HPP

#include "siconos_storage.hpp"
#include "siconos_pattern.hpp"
#include "siconos_linear_algebra.hpp"
#include <range/v3/view/zip.hpp>

#define GET(X) decltype(auto) X() \
  { return get<typename Handle::type::X>(*default_interface<Handle>::self(), \
                                         default_interface<Handle>::self()->data());}

#define ITEM(Y)                                         \
  decltype(auto) Y()                                    \
  {                                                     \
    return make_full_handle<Y##_t>(self()->get(), self()->data()); \
  };
namespace siconos
{

  decltype(auto) edge1(auto& g, auto& descr)
  {
    auto [oei, oee] = g.out_edges(descr);
    return oei;
  };

  decltype(auto) edge2(auto& g, auto& descr)
  {
    auto [oei, oee] = g.out_edges(descr);
    return oee;
  };

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

      using attributes = type::attributes<mass_matrix, q, velocity, fext>;

      template<typename Handle>
      struct interface : default_interface<Handle>
      {
        using default_interface<Handle>::self;

        decltype(auto) mass_matrix()
        { return Handle::type::mass_matrix::at(*self()); };

        GET(velocity);
        GET(q);
        GET(fext);
      };
    };

    struct relation : item<>
    {
      struct h_matrix : some::matrix<dof, dof>, access<h_matrix> {};

      using attributes = type::attributes<h_matrix>;

      using properties = std::tuple<attached_storage<dynamical_system, symbol<"p0">,
                                                     some::vector<dof>>>;

      template<typename Handle>
      struct interface : default_interface<Handle>
      {

        using default_interface<Handle>::self;

        GET(h_matrix);

        template<typename Interaction, typename DS>
        void compute_input
        (double time, Interaction inter, DS ds, auto level)
        {
          auto& lambda = inter.lambda()[level];
          prod(lambda, self()->h_matrix(), ds.property(symbol<"p0">{}));
         // link to ds variables => check graph
        };

        template<typename Interaction, typename DS>
        void compute_output(double time, Interaction inter, DS ds, auto level)
        {
          auto& y = inter.y()[1];
          auto& H = self()->h_matrix();

          auto& x = [&ds,&level]() -> decltype(auto)
          {
          if constexpr (level == 1)
          {
            return ds.velocity();
          }
          else
          {
            return ds.q();
          }
          }();

          prod(H, x, y);
        };
      };
    };
  };

  struct equality_condition_nsl {};
  struct relay_nsl {};
  struct nonsmooth_law
  {
    struct newton_impact_friction : item<>
    {
      struct e : access<e>, some::scalar {};
      struct mu : some::scalar, access<mu> {};

      using attributes = type::attributes<e, mu>;

      template<typename Handle>
      struct interface : default_interface<Handle>
      {
        using default_interface<Handle>::self;

        GET(e);
        GET(mu);

      };

    };

    struct newton_impact : item<>
    {
      struct e : some::scalar, access<e>, text<"e"> {};

      using attributes = gather<e>;

      template<typename Handle>
      struct interface : default_interface<Handle>
      {
        using default_interface<Handle>::self;
        GET(e);
      };
    };
  };

  template<typename Formulation, typename Interaction>
  struct topology : item<>
  {
    using dynamical_system = typename Formulation::dynamical_system;
    using interaction = Interaction;

    struct dynamical_system_graphs :
      some::bounded_collection<some::graph<some::item_ref<dynamical_system>,
                                           some::item_ref<interaction>>, 1>{};

    struct interaction_graphs :
      some::bounded_collection<some::graph<some::item_ref<interaction>,
                                           some::item_ref<dynamical_system>>, 2>{};


    using attributes = gather<dynamical_system_graphs, interaction_graphs>;

    template<typename Handle>
    struct interface : default_interface<Handle>
    {
      GET(dynamical_system_graphs);
      GET(interaction_graphs);

      using default_interface<Handle>::self;
      template <match::handle<dynamical_system> ds_t,
        match::handle<interaction> inter_t>
      decltype(auto) link (inter_t inter, ds_t ds)
      {
        auto& dsg0 = self()->dynamical_system_graphs()[0];
        auto& ig0 = self()->interaction_graphs()[0];;

        auto dsgv = dsg0.add_vertex(ds);
        return dsg0.add_edge(dsgv, dsgv, inter, ig0);
      };
    };
  };

  template<typename Nslaw, typename Formulation, std::size_t N, std::size_t K=2>
  struct interaction : item<>
  {
    struct nonsmooth_law : some::item_ref<Nslaw>, access<nonsmooth_law> {};
    struct relation      : some::item_ref<typename Formulation::relation>,
                           access<relation> {};
    static constexpr auto number_of_ds = N;
    struct lambda : some::vector<K, some::vector<N*Formulation::dof>>, access<lambda> {};
    struct y : some::vector<K, some::vector<N*Formulation::dof>>, access<y> {};
    using attributes = gather<nonsmooth_law,
                              relation, lambda, y>;

    template<typename Handle>
    struct interface : default_interface<Handle>
    {
      using default_interface<Handle>::self;

      GET(nonsmooth_law);
      GET(lambda);
      auto relation()
      { return full_handle(Handle::type::relation::at(*self()), self()->data()); };
      GET(y);
    };

  };

  struct lcp
  {};

  template<typename Type>
  struct one_step_nonsmooth_problem : item<>
  {
    using problem_type = Type;
    struct level : some::indice {};
    using attributes = gather<level>;

    template<typename Handle>
    struct interface : default_interface<Handle>
    {
      using default_interface<Handle>::self;

      GET(level);
    };

  };

  template<typename Form, typename Interaction>
  struct one_step_integrator
  {
    using formulation = Form;
    using interaction = Interaction;
    using system = typename formulation::dynamical_system;
    using q = typename system::q;
    using velocity = typename system::velocity;
    using mass_matrix = typename system::mass_matrix;
    using fext = typename system::fext;

    struct euler : item<>
    {
      using properties = gather<
        keep<q, 2>,
        keep<velocity, 2>>;

      using attributes = gather<>;

      template<typename Handle>
      struct interface : default_interface<Handle>
      {
        using default_interface<Handle>::self;
        void compute_free_state(auto step, auto h)
        {
          // auto& data = self()->data();
          // auto& velocities = get_memory<velocity>(data);
          // auto& mass_matrices = get_memory<mass_matrix>(data);
          // auto& external_forces = get_memory<fext>(data);

          // auto& Ms =      memory(step, mass_matrices);
          // auto& vs =      memory(step, velocities);
          // auto& vs_next = memory(step+1, velocities);
          // auto& fs =      memory(step, external_forces);

        };
      };

    };

    struct moreau_jean : item<>
    {
      struct theta : some::scalar {};
      struct gamma : some::scalar {};
      struct constraint_activation_threshold : some::scalar {};
      using attributes = type::attributes<theta, gamma,
                                          constraint_activation_threshold>;

      using properties = gather<
        keep<typename system::q, 2>,
        keep<typename system::velocity, 2>>;

      template<typename Handle>
      struct interface : default_interface<Handle>
      {
        using default_interface<Handle>::self;

        GET(theta);
        GET(gamma);
        GET(constraint_activation_threshold);

        bool add_interaction_in_index_set(
          auto inter,
          auto h,
          auto i)
        {
          auto          y = inter.y()[i-1][0];
          const auto ydot = inter.y()[i][0];
          const auto& gamma_v = 0.5;
          y += gamma_v * h * ydot ;

          return y <= self()->constraint_activation_threshold();
        };


        bool remove_interaction_from_index_set(
          auto inter,
          auto h,
          auto i)
        {
          return !add_interaction_in_index_set(inter, h, i);
        };

        auto compute_iteration_matrix(auto step)
        {
          auto& data = self()->data();
          auto& mass_matrices = get_memory<mass_matrix>(data);
          auto& external_forces = get_memory<fext>(data);

          auto& mats =      memory(step, mass_matrices);
          auto& fs =      memory(step, external_forces);

          if constexpr(has_property<mass_matrix, some::time_invariant>(data))
          {
            if constexpr(has_property<fext, some::time_invariant>(data))
            {
              if constexpr(has_property<mass_matrix, some::diagonal>(data))
              {
                for (auto [mat, f] : ranges::views::zip(mats, fs))
                {
                  linear_algebra::solve_in_place(mat, f);
                }
              }
            }
          }
        };
        auto compute_free_state(auto step, auto h)
        {
          auto& data = self()->data();
          auto& velocities = get_memory<velocity>(data);
          auto fexts  = get_memory<fext>(data);
          auto theta_ = self()->theta();

          auto& vs =      memory(step, velocities);
          auto& vs_next = memory(step+1, velocities);
          auto& minv_fs      = memory(step, fexts);
          auto& minv_fs_next      = memory(step+1, fexts);

          // f <- M  f
          for (auto [v, v_next, minv_f, minv_f_next] :
                 ranges::views::zip(vs, vs_next, minv_fs, minv_fs_next))
          {
            v_next = v + h*theta_*minv_f + h*(1 - theta_)*minv_f_next;
          }
        };
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

    template<typename Handle>
    struct interface : default_interface<Handle>
    {
      using default_interface<Handle>::self;

      GET(t0);
      GET(h);
      GET(step);
    };

  };

  template<typename TD, typename OSI, typename OSNSPB, typename TOPO>
  struct time_stepping : item<>
  {
    using time_discretization_t = TD;
    using one_step_integrator_t = OSI;
    using one_step_nonsmooth_problem_t = OSNSPB;
    using topology_t = TOPO;

    using attributes = type::attributes_of_items<time_discretization_t,
                                                 one_step_integrator_t,
                                                 one_step_nonsmooth_problem_t,
                                                 topology_t>;

    template<typename Handle>
    struct interface : default_interface<Handle>
    {
      using default_interface<Handle>::self;

      ITEM(time_discretization);
      ITEM(one_step_integrator);
      ITEM(one_step_nonsmooth_problem);
      ITEM(topology);

      decltype(auto) current_step()
      {
        return (*self()).time_discretization().step();
      }

      decltype(auto) time_step()
      {
        return self()->time_discretization().h();
      }

      decltype(auto) compute_output(auto level)
      {
        auto& index_set1 = topology().interaction_graphs()[1];

        for (auto [ui1, ui1end] = index_set1.vertices();
             ui1 != ui1end; ++ui1)
        {
          auto inter1 = full_handle(index_set1.bundle(*ui1), self()->data());
          auto [oei, oee] = index_set1.out_edges(*ui1);
          if (oee == oei)
          {
            // one ds
            inter1.relation().compute_output(0.,
                                             inter1,
                                             full_handle(index_set1.bundle(*oei), self()->data()),
                                             level);
          };
        }
      }
      void compute_one_step()
      {
        self()->one_step_integrator().compute_free_state(current_step(),
                                                         time_step());
        current_step() += 1;
      }

      void update_indexsets(auto i)
      {
        auto& data = self()->data();
        using info_t = std::decay_t<decltype(ground::get<info>(data))>;
        using env = typename info_t::env;

        auto osi = self()->one_step_integrator();
        auto topo = self()->topology();

        auto& index_set0 = topo.interaction_graphs()[0];
        auto& index_set1 = topo.interaction_graphs()[1];
//        auto& dsg0 = topo.dynamical_system_graphs()[0];

        // Check index_set1
        auto [ui1, ui1end] = index_set1.vertices();

        // Remove interactions from the index_set1
        for (auto v1next = ui1; ui1 != ui1end; ui1 = v1next)
        {
          ++v1next;
          auto inter1 = full_handle(index_set1.bundle(*ui1), self()->data());  // get inter handle
//          auto rel1 = inter1.relation();

          if (index_set0.is_vertex(inter1))
          {
            auto inter1_descr0 = index_set0.descriptor(inter1);
            assert((index_set0.color(inter1_descr0) == env::white_color));

            index_set0.color(inter1_descr0) = env::gray_color;
            if constexpr (
              !std::derived_from<typename topology_t::interaction::nonsmooth_law,
              equality_condition_nsl>)
            {
              // We assume that the integrator of the ds1 drive the update of the index set
              //SP::OneStepIntegrator Osi = index_set1.properties(*ui1).osi;
//              auto&& ds1 = edge1(index_set1, *ui1);
//              auto& osi = dsg0.properties(dsg0.descriptor(ds1)).osi;

//              //if(predictorDeactivate(inter1,i))
              if(osi.remove_interaction_from_index_set(inter1, self()->time_step(), i))
              {
//                // Interaction is not active
//                // ui1 becomes invalid
                index_set0.color(inter1_descr0) = env::black_color;
//                index_set1.eraseProperties(*ui1);

//              InteractionsGraph::OEIterator oei, oeiend;
                for(auto [oei, oeiend] = index_set1.out_edges(*ui1);
                    oei != oeiend; ++oei)
                {
                  auto [ed1, ed2] = index_set1.edges(index_set1.source(*oei), index_set1.target(*oei));
                  if(ed2 != ed1)
                  {
//                    index_set1.eraseProperties(ed1);
//                    index_set1.eraseProperties(ed2);
                  }
                  else
                  {
//                    index_set1.eraseProperties(ed1);
                  }
                }
                index_set1.remove_vertex(inter1);
//           /* \warning V.A. 25/05/2012 : Multiplier lambda are only set to zero if they are removed from the IndexSet*/
                inter1.lambda()[1] = {};
//                topo->setHasChanged(true);
              }
            }
          }
          else
          {
            // Interaction is not in index_set0 anymore.
            // ui1 becomes invalid
//            index_set1.eraseProperties(*ui1);
            for(auto [oei, oeiend] = index_set1.out_edges(*ui1);
                oei != oeiend; ++oei)
            {
              auto [ed1, ed2] = index_set1.edges(index_set1.source(*oei), index_set1.target(*oei));
              if(ed2 != ed1)
              {
//                index_set1.eraseProperties(ed1);
//                index_set1.eraseProperties(ed2);
              }
              else
              {
//                index_set1.eraseProperties(ed1);
              }
            }

            index_set1.remove_vertex(inter1);
//       topo->setHasChanged(true);
          }
      }

//   // index_set0\index_set1 scan
//   InteractionsGraph::VIterator ui0, ui0end;
//   //Add interaction in index_set1
      for(auto [ui0, ui0end] = index_set0.vertices(); ui0 != ui0end; ++ui0)
      {
        if(index_set0.color(*ui0) == env::black_color)
        {
          // reset
          index_set0.color(*ui0) = env::white_color ;
        }
        else
        {
          if(index_set0.color(*ui0) == env::gray_color)
          {
            // reset
            index_set0.color(*ui0) = env::white_color;

            assert(index_set1.is_vertex(index_set0.bundle(*ui0)));
//         /*assert( { !predictorDeactivate(index_set0->bundle(*ui0),i) ||
//           Type::value(*(index_set0->bundle(*ui0)->nonSmoothLaw())) == Type::EqualityConditionNSL ;
//           });*/
          }
          else
          {
            assert(index_set0.color(*ui0) == env::white_color);

            auto inter0 = full_handle(index_set0.bundle(*ui0), self()->data());
            assert(!index_set1.is_vertex(inter0));
            bool activate = true;
            if constexpr(
              !std::derived_from<typename topology_t::interaction::nonsmooth_law,
              equality_condition_nsl> &&
              !std::derived_from<typename topology_t::interaction::nonsmooth_law,
              relay_nsl>)
//             && Type::value(*(inter0->nonSmoothLaw())) != Type::RelayNSL)
            {
              //SP::OneStepIntegrator Osi = index_set0->properties(*ui0).osi;
//           // We assume that the integrator of the ds1 drive the update of the index set
//              auto&& ds1 = edge1(index_set1, *ui0);
//           OneStepIntegrator& osi = *DSG0.properties(DSG0.descriptor(ds1)).osi;

              activate = osi.add_interaction_in_index_set(inter0, time_step(), i);
            }
            if(activate)
            {
              assert(!index_set1.is_vertex(inter0));

//           // vertex and edges insertion in index_set1
              index_set1.copy_vertex(inter0, index_set0);
//           topo->setHasChanged(true);
              assert(index_set1.is_vertex(inter0));
            }
          }
        }
      }

//   assert(index_set1->size() <= index_set0->size());

//   DEBUG_PRINTF("TimeStepping::updateIndexSet(unsigned int i). update index_sets end : index_set0 size : %ld\n", index_set0->size());
//   DEBUG_PRINTF("TimeStepping::updateIndexSet(unsigned int i). update IndexSets end : index_set1 size : %ld\n", index_set1->size());
// }

      };
    };
  };
};
#endif
