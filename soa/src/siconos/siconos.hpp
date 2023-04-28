#pragma once

#include "siconos/storage/storage.hpp"
#include "siconos/utils/pattern.hpp"
#include "siconos/algebra/linear_algebra.hpp"
#include "siconos/algebra/numerics.hpp"
#include <range/v3/view/zip.hpp>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/filter.hpp>
#include <tuple>

#define GET(X) decltype(auto) X()                                       \
  { return get<typename Handle::type::X>(*default_interface<Handle>::self(), \
                                         default_interface<Handle>::self()->data());}

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
        some::matrix<some::scalar, dof, dof>,
        access<mass_matrix> {};

      struct q : some::vector<some::scalar, dof>,
                 access<q> {};

      struct velocity : some::vector<some::scalar, dof>,
                        access<velocity> {};

      // should not be an attribute
      struct fext : some::vector<some::scalar, dof>, // some::function<...>
                    access<fext> {};

      using attributes = types::attributes<mass_matrix, q, velocity, fext>;

      template<typename Handle>
      struct interface : default_interface<Handle>
      {
        using default_interface<Handle>::self;

        decltype(auto) mass_matrix()
        { return Handle::type::mass_matrix::at(*self()); }

        decltype(auto) velocity()
        { return Handle::type::velocity::at(*self()); }

        decltype(auto) q()
        { return Handle::type::q::at(*self()); }

        decltype(auto) fext()
        { return Handle::type::fext::at(*self()); }
      };
    };

    struct relation : item<>
    {

      using attributes = types::attributes<>;

      template<typename Handle>
      struct interface : default_interface<Handle>
      {

        using default_interface<Handle>::self;

        template<match::any_full_handle DynamicalSystem,
                 match::any_full_handle Interaction>
        void compute_input
        (double time, DynamicalSystem ds, Interaction inter, auto level)
        {
          auto& lambda = inter.lambda()[level][0];
          prod(trans(inter.h_matrix()[0]), lambda, ds.property(symbol<"p0">{})[0]);
         // link to ds variables => check graph
        };

        template<match::any_full_handle DynamicalSystem,
                 match::any_full_handle Interaction>
        void compute_output(double time, DynamicalSystem ds,
                            Interaction inter, auto level)
        {
          auto& y = inter.y()[1][0];
          auto& H = inter.h_matrix()[0];

          auto& x = [&level,&ds]() -> decltype(auto)
          {
          if constexpr (level == 1)
          {
            return ds.velocity();
          }
          else
          {
            return ds.q();;
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

      static constexpr auto size = 2;
      using attributes = types::attributes<e, mu>;

      template<typename Handle>
      struct interface : default_interface<Handle>
      {
        using default_interface<Handle>::self;


        GET(e)
        GET(mu)

      };

    };

    struct newton_impact : item<>
    {
      struct e : some::scalar, access<e>, text<"e"> {};
      static constexpr auto size = 1;
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

    // number of involved ds
    struct nids : some::indice {};

    using attributes = gather<
      dynamical_system_graphs,
      interaction_graphs, nids>;

    using properties = gather<
      attached_storage<dynamical_system, symbol<"involved">, some::boolean>,
        attached_storage<dynamical_system, symbol<"index">, some::indice>,
        attached_storage<
            dynamical_system, symbol<"p0">,
            some::vector<some::vector<some::scalar, Formulation::dof>, 2>>,
        attached_storage<
            dynamical_system, symbol<"velocity">,
            some::vector<some::vector<some::scalar, Formulation::dof>, 2>>,
        attached_storage<
            dynamical_system, symbol<"q">,
            some::vector<some::vector<some::scalar, Formulation::dof>, 2>>,
        attached_storage<interaction, symbol<"nds">, some::indice>,
        attached_storage<interaction, symbol<"ds1">,
                         some::item_ref<dynamical_system>>,
        attached_storage<interaction, symbol<"ds2">,
                         some::item_ref<dynamical_system>>,
        attached_storage<
            dynamical_system, symbol<"vd">,
            some::vdescriptor<typename dynamical_system_graphs::type>>,
        attached_storage<
            interaction, symbol<"vd">,
            some::vdescriptor<typename interaction_graphs::type>>>;

    template<typename Handle>
    struct interface : default_interface<Handle>
    {
      GET(dynamical_system_graphs);
      GET(interaction_graphs);
      GET(nids);

      using default_interface<Handle>::self;

      template<match::handle<dynamical_system> Hds>
      decltype(auto) link (Hds ds)
      {
        auto& data = self()->data();
        auto& dsg0 = self()->dynamical_system_graphs()[0];
        auto& ig0 = self()->interaction_graphs()[0];;

        auto inter = add<interaction>(data);
        auto dsgv = dsg0.add_vertex(ds);
        auto [dsg0_ed, ig0_vd] = dsg0.add_edge(dsgv, dsgv, inter, ig0);
        ds.property(symbol<"vd">{}) = dsgv;
        inter.property(symbol<"vd">{}) = ig0_vd;
        inter.property(symbol<"nds">{}) = 1;
        inter.template property<"ds1">() = ds;

        return inter;
      };

      template<match::handle<dynamical_system> Hds>
      decltype(auto) link (Hds ds1, Hds ds2)
      {
        auto& data = self()->data();
        auto& dsg0 = self()->dynamical_system_graphs()[0];
        auto& ig0 = self()->interaction_graphs()[0];;

        auto inter = add<interaction>(data);
        auto dsgv1 = dsg0.add_vertex(ds1);
        auto dsgv2 = dsg0.add_vertex(ds2);
        auto [dsg0_ed, ig0_vd] = dsg0.add_edge(dsgv1, dsgv2, inter, ig0);

        ds1.property(symbol<"vd">{}) = dsgv1;
        ds2.property(symbol<"vd">{}) = dsgv2;

        inter.property(symbol<"nds">{}) = 2;
        inter.property(symbol<"vd">{}) = ig0_vd;
        inter.property(symbol<"ds1">{}) = ds1;
        inter.property(symbol<"ds2">{}) = ds2;

        return inter;
      };

      // strategy 2 : h_matrix is assembled for involved ds only
      auto make_index()
      {
        auto& data = self()->data();
        using info_t = std::decay_t<decltype(ground::get<info>(data))>;
        using env = typename info_t::env;
        using indice = typename env::indice;
        auto& dsg0 = self()->dynamical_system_graphs()[0];
        indice counter = 0;
        for (auto [dsi, dsiend] = dsg0.vertices(); dsi!=dsiend; ++dsi)
        {
          indice& prop = handle(dsg0.bundle(*dsi), data).property(symbol<"index">{});
          auto&& involved = handle(dsg0.bundle(*dsi), data).property(symbol<"involved">{});
          auto [oei, oeiend] = dsg0.out_edges(*dsi);
          if (oeiend != oei) // ds is involved in some interaction
          {
            involved = true;
            prop = counter++;
          } else
          {
            involved = false;
            prop = 0;
          }
        }
        self()->nids() = counter;
      };
    };
  };

  template<typename Nslaw, typename Formulation, std::size_t K=2>
  struct interaction : item<>
  {
    struct nonsmooth_law : some::item_ref<Nslaw>, access<nonsmooth_law> {};
    struct relation      : some::item_ref<typename Formulation::relation>,
                           access<relation> {};

    struct h_matrix
        : some::vector<
              some::matrix<some::scalar, Formulation::dof, Nslaw::size>, 2>,
          access<h_matrix> {};
    struct lambda
      : some::vector<some::vector<some::vector<some::scalar, Nslaw::size>, 2>, 2>,
          access<lambda> {};
    struct y
      : some::vector<some::vector<some::vector<some::scalar, Nslaw::size>, 2>, 2>,
          access<y> {};

    using attributes = gather<nonsmooth_law, relation, h_matrix, lambda, y>;

    template<typename Handle>
    struct interface : default_interface<Handle>
    {
      using default_interface<Handle>::self;

      GET(nonsmooth_law)
      GET(lambda)
      decltype(auto) relation()
      { return handle(Handle::type::relation::at(*self()), self()->data()); }

      decltype(auto) h_matrix()
      { return Handle::type::h_matrix::at(*self()); }

      GET(y)
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
    using y = typename interaction::y;
    using q = typename system::q;
    using velocity = typename system::velocity;
    using mass_matrix = typename system::mass_matrix;
    using h_matrix = typename interaction::h_matrix;
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

      struct h_matrix_assembled : some::unbounded_matrix<h_matrix> {};
      struct q_vector_assembled : some::unbounded_vector<q> {};
      struct velocity_vector_assembled : some::unbounded_matrix<velocity> {};
      struct y_vector_assembled : some::unbounded_matrix<y> {};
      struct mass_matrix_assembled : some::unbounded_matrix<mass_matrix> {};
      struct w_matrix : some::unbounded_matrix<mass_matrix> {};

      using attributes = types::attributes<theta,
                                           gamma,
                                           constraint_activation_threshold,
                                           h_matrix_assembled,
                                           mass_matrix_assembled,
                                           w_matrix,
                                           q_vector_assembled,
                                           velocity_vector_assembled,
                                           y_vector_assembled>;

      using properties = gather<
        keep<typename system::q, 2>,
        keep<typename system::velocity, 2>>;

      template<typename Handle>
      struct interface : default_interface<Handle>
      {
        using default_interface<Handle>::self;

        GET(theta)
        GET(gamma)
        GET(constraint_activation_threshold)
        GET(h_matrix_assembled)
        GET(q_vector_assembled)
        GET(velocity_vector_assembled)
        GET(y_vector_assembled)
        GET(mass_matrix_assembled)
        GET(w_matrix)

        bool add_interaction_in_index_set(
          auto inter,
          auto h,
          auto i)
        {
          auto          y = inter.y()[i-1][0][0];
          const auto ydot = inter.y()[i][0][0];
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
        }

        // strategy 1 : assemble the matrix for involved ds only
        auto assemble_h_matrix_for_involved_ds(auto step)
        {
          auto& data = self()->data();
          auto& h_matrices = memory(step, get_memory<h_matrix>(data));
          auto size = std::size(h_matrices);
          auto& ids1s =
              ground::get<attached_storage<interaction, symbol<"ds1">,
                                           some::item_ref<system>>>(data)[0];
          auto& ids2s =
              ground::get<attached_storage<interaction, symbol<"ds2">,
                                           some::item_ref<system>>>(data)[0];

          self()->h_matrix_assembled().resize(size, size, false);
          for (auto [mat, ids1, ids2] :
               ranges::views::zip(h_matrices, ids1s, ids2s)) {
            auto i = handle(ids1, data)
                         .property(symbol<"index">{});  // involved index
            auto j = handle(ids2, data)
                         .property(symbol<"index">{});  // involved index
            self()->h_matrix_assembled()(i, j) = mat;
          }
        }

        // strategy 2 : assemble the whole matrix (size = number of ds)
        auto assemble_h_matrix_for_all_ds(auto step)
        {
          auto& data = self()->data();

          // size is the number of ds
          auto size = std::size(memory(step, get_memory<mass_matrix>(data)));

          auto& h_matrices = memory(step, get_memory<h_matrix>(data));
          auto& ids1s =
              ground::get<attached_storage<interaction, symbol<"ds1">,
                                           some::item_ref<system>>>(data)[0];
          auto& ids2s =
              ground::get<attached_storage<interaction, symbol<"ds2">,
                                           some::item_ref<system>>>(data)[0];

          self()->h_matrix_assembled().resize(size, size, false);
          for (auto [mat, ids1, ids2] :
               ranges::views::zip(h_matrices, ids1s, ids2s)) {
            auto i = handle(ids1, data).get();  // global index
            auto j = handle(ids2, data).get();  // global index
            self()->h_matrix_assembled()(i, j) = mat;
          }
        }

        auto assemble_mass_matrix_for_involved_ds(auto step)
        {
          auto& data = self()->data();

          auto size = std::size(memory(step, get_memory<h_matrix>(data)));

          auto& mass_matrices = memory(step, get_memory<mass_matrix>(data));
          auto& involved_ds = ground::get<attached_storage<interaction, symbol<"involved">,
                                                         some::item_ref<system>>>(data)[0];

          self()->mass_matrix_assembled().resize(size, size);

          for (auto [i, mat] : ranges::views::enumerate(mass_matrices) |
                                   ranges::views::filter([](auto k_m) {
                                     auto [k,_] = k_m;
                                     return involved_ds[k];
                                   }))
          {
            self()->mass_matrix_assembled()(i, i) = mat;
          }
        }

        auto assemble_mass_matrix_for_all_ds(auto step)
        {
          auto& data = self()->data();

          // size is the number of ds
          auto size = std::size(memory(step, get_memory<mass_matrix>(data)));

          auto& mass_matrices = memory(step, get_memory<mass_matrix>(data));

          self()->mass_matrix_assembled().resize(size, size);

          // also mapping block vector -> block diagonal matrix
          for (auto [i, mat] : ranges::views::enumerate(mass_matrices)) {
            self()->mass_matrix_assembled()(i,i) = mat;
          }
        }

        auto compute_w_matrix(auto step)
        {
          auto tmp_matrix = std::decay_t<decltype(w_matrix())> {};
          solve(mass_matrix_assembled(),
                h_matrix_assembled(), tmp_matrix);
          // aliasing ?
          prod(trans(h_matrix_assembled()), tmp_matrix, w_matrix());
        }

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
            // beware of temporaries...
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

  template<match::item ...Items>
  struct time_stepping : item<>
  {
    using items = gather<Items...>;
    using time_discretization_t = types::nth_t<0, items>;
    using one_step_integrator_t = types::nth_t<1, items>;
    using one_step_nonsmooth_problem_t = types::nth_t<2, items>;
    using topology_t = types::nth_t<3, items>;

    using attributes = types::attributes_of_items<Items...>;

    template<typename Handle>
    struct interface : default_interface<Handle>
    {
      using default_interface<Handle>::self;

      decltype(auto) time_discretization()
      {
        return make_handle<time_discretization_t>(*self());
      }
      decltype(auto) one_step_integrator()
      {
        return make_handle<one_step_integrator_t>(*self());
      }
      decltype(auto) one_step_nonsmooth_problem()
      {
        return make_handle<one_step_nonsmooth_problem_t>(*self());
      }
      decltype(auto) topology()
      {
        return make_handle<topology_t>(*self());
      }
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
        auto osi = one_step_integrator();

        switch (level)
        {
        case 1 :
          prod(osi.h_matrix_assembled(),
               osi.velocity_vector_assembled(),
               osi.y_vector_assembled());
          break;
        case 0 :
          prod(osi.h_matrix_assembled(),
               osi.q_vector_assembled(),
               osi.y_vector_assembled());
          break;
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
          auto inter1 = handle(index_set1.bundle(*ui1), self()->data());  // get inter handle
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

            auto inter0 = handle(index_set0.bundle(*ui0), self()->data());
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

