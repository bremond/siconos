#pragma once

#include "siconos/algebra/eigen.hpp"
#include "siconos/model/nslaws.hpp"
#include "siconos/simul/simul.hpp"

namespace siconos::simul {
template <match::item... Items>
struct time_stepping : item<> {
  using items = gather<Items...>;
  using time_discretization_t = nth_t<0, items>;
  using one_step_integrator_t = nth_t<1, items>;
  using one_step_nonsmooth_problem_t = nth_t<2, items>;
  using topology_t = nth_t<3, items>;

  using formulation_t =
      typename one_step_nonsmooth_problem_t::problem_t::formulation_t;
  using attributes = types::attributes_of_items<Items...>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
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
    decltype(auto) topology() { return make_handle<topology_t>(*self()); }
    decltype(auto) current_step()
    {
      return (*self()).time_discretization().step();
    }

    decltype(auto) time_step() { return self()->time_discretization().h(); }

    auto compute_one_step()
    {
      auto osi = one_step_integrator();
      auto step = current_step();

      osi.compute_free_state(step, time_step());

      // -> ydot (step+1)
      osi.update_positions(step, time_step());
      osi.compute_output(step);
      osi.compute_output(step + 1);

      // activations of interactions
      update_indexsets(0);

      // compute active interactions
      auto [ninter, nds] = osi.compute_active_interactions(step, time_step());

      if (nds > 0) {
        // a least one activated interaction

        //        print("ninter, nds = {},{}\n", ninter, nds);

        osi.assemble_h_matrix_for_involved_ds(step, ninter, nds);
        osi.assemble_mass_matrix_for_involved_ds(step, nds);

        osi.resize_assembled_vectors(step, ninter);

        // H M^-1 H^t
        osi.compute_w_matrix(step);
        osi.nsl_effect_on_free_output(step);
        osi.compute_q_nsp_vector_assembled(step, ninter);

        self()->template solve_nonsmooth_problem<formulation_t>(step, ninter);

        osi.compute_input();
        osi.update_all_velocities(step);
        osi.update_positions(step, time_step());
      }
      else {
        print(".");
      }

      // do nothing if lagrangian_r is time_invariant
      // osi.update_h_matrices(step);

      // do nothing if fext is time_invariant
      // osi.update_iteration_matrix(step);

      current_step() += 1;

      print("step {}\n", current_step());
      return nds;
    }

    template <typename Formulation>
    void solve_nonsmooth_problem(auto step, auto ninter)
    {  // Mz=w+q
      auto osi = self()->one_step_integrator();

      //      resize(osi.lambda_vector_assembled(), ninter);
      //      resize(osi.ydot_vector_assembled(), ninter);

      self()->one_step_nonsmooth_problem().template solve<Formulation>(
          osi.w_matrix(),                 // M
          osi.q_nsp_vector_assembled(),   // q
          osi.lambda_vector_assembled(),  // z
          osi.ydot_vector_assembled());   // w
    }

    bool has_next_event()
    {
      return time_discretization().t0() +
                 current_step() * time_discretization().h() <
             time_discretization().tmax();
    }
    void update_indexsets(auto i)
    {
      auto& data = self()->data();
      using info_t = std::decay_t<decltype(ground::get<storage::info>(data))>;
      using env = typename info_t::env;

      auto topo = self()->topology();

      auto& index_set0 = topo.interaction_graphs()[0];
      auto& index_set1 = topo.interaction_graphs()[1];
      //        auto& dsg0 = topo.dynamical_system_graphs()[0];

      // Check index_set1
      auto [ui1, ui1end] = index_set1.vertices();

      // Remove interactions from the index_set1
      for (auto v1next = ui1; ui1 != ui1end; ui1 = v1next) {
        ++v1next;
        auto inter1 = storage::handle(
            self()->data(), index_set1.bundle(*ui1));  // get inter handle
        //          auto rel1 = inter1.relation();

        if (index_set0.is_vertex(inter1)) {
          auto inter1_descr0 = index_set0.descriptor(inter1);
          assert((index_set0.color(inter1_descr0) == env::white_color));

          index_set0.color(inter1_descr0) = env::gray_color;
          if constexpr (!std::derived_from<
                            attr_t<typename topology_t::interaction, "nslaw">,
                            model::equality_condition>) {
            // We assume that the integrator of the ds1 drive the update of
            // the index set SP::OneStepIntegrator Osi =
            // index_set1.properties(*ui1).osi;
            //              auto&& ds1 = edge1(index_set1, *ui1);
            //              auto& osi =
            //              dsg0.properties(dsg0.descriptor(ds1)).osi;

            //              //if(predictorDeactivate(inter1,i))
            //            if (osi.remove_interaction_from_index_set(
            if (!inter1.property(symbol<"activation">{})) {
              //                // Interaction is not active
              //                // ui1 becomes invalid
              index_set0.color(inter1_descr0) = env::black_color;
              //                index_set1.eraseProperties(*ui1);

              //              InteractionsGraph::OEIterator oei, oeiend;
              for (auto [oei, oeiend] = index_set1.out_edges(*ui1);
                   oei != oeiend; ++oei) {
                auto [ed1, ed2] = index_set1.edges(index_set1.source(*oei),
                                                   index_set1.target(*oei));
                if (ed2 != ed1) {
                  //                    index_set1.eraseProperties(ed1);
                  //                    index_set1.eraseProperties(ed2);
                }
                else {
                  //                    index_set1.eraseProperties(ed1);
                }
              }
              index_set1.remove_vertex(inter1);
              //           /* \warning V.A. 25/05/2012 : Multiplier lambda are
              //           only set to zero if they are removed from the
              //           IndexSet*/
              inter1.lambda() = {};
              //                topo->setHasChanged(true);
            }
          }
        }
        else {
          // Interaction is not in index_set0 anymore.
          // ui1 becomes invalid
          //            index_set1.eraseProperties(*ui1);
          for (auto [oei, oeiend] = index_set1.out_edges(*ui1); oei != oeiend;
               ++oei) {
            auto [ed1, ed2] = index_set1.edges(index_set1.source(*oei),
                                               index_set1.target(*oei));
            if (ed2 != ed1) {
              //                index_set1.eraseProperties(ed1);
              //                index_set1.eraseProperties(ed2);
            }
            else {
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
      for (auto [ui0, ui0end] = index_set0.vertices(); ui0 != ui0end; ++ui0) {
        if (index_set0.color(*ui0) == env::black_color) {
          // reset
          index_set0.color(*ui0) = env::white_color;
        }
        else {
          if (index_set0.color(*ui0) == env::gray_color) {
            // reset
            index_set0.color(*ui0) = env::white_color;

            assert(index_set1.is_vertex(index_set0.bundle(*ui0)));
            //         /*assert( {
            //         !predictorDeactivate(index_set0->bundle(*ui0),i) ||
            //           Type::value(*(index_set0->bundle(*ui0)->nonSmoothLaw()))
            //           == Type::EqualityConditionNSL ;
            //           });*/
          }
          else {
            assert(index_set0.color(*ui0) == env::white_color);

            auto inter0 =
                storage::handle(self()->data(), index_set0.bundle(*ui0));
            assert(!index_set1.is_vertex(inter0));
            bool activate = true;
            if constexpr (
                !std::derived_from<
                    attr_t<typename topology_t::interaction, "nslaw">,
                    model::equality_condition> &&
                !std::derived_from<
                    attr_t<typename topology_t::interaction, "nslaw">,
                    model::relay>)
            //             && Type::value(*(inter0->nonSmoothLaw())) !=
            //             Type::RelayNSL)
            {
              // SP::OneStepIntegrator Osi = index_set0->properties(*ui0).osi;
              //           // We assume that the integrator of the ds1 drive
              //           the update of the index set
              //              auto&& ds1 = edge1(index_set1, *ui0);
              //           OneStepIntegrator& osi =
              //           *DSG0.properties(DSG0.descriptor(ds1)).osi;
              activate = inter0.property(symbol<"activation">{});
              //                  osi.add_interaction_in_index_set(inter0,
              //                  time_step(), i);
            }
            if (activate) {
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

      //   DEBUG_PRINTF("TimeStepping::updateIndexSet(unsigned int i). update
      //   index_sets end : index_set0 size : %ld\n", index_set0->size());
      //   DEBUG_PRINTF("TimeStepping::updateIndexSet(unsigned int i). update
      //   IndexSets end : index_set1 size : %ld\n", index_set1->size());
      // }
    };

    void initialize() { one_step_integrator().initialize(current_step()); }

    auto methods()
    {
      using env_t = decltype(self()->env());
      using indice = typename env_t::indice;
      //      using scalar = typename env_t::scalar;

      return collect(
          method("time_discretization",
                 &interface<Handle>::time_discretization),
          method("one_step_integrator",
                 &interface<Handle>::one_step_integrator),
          method("one_step_nonsmooth_problem",
                 &interface<Handle>::one_step_nonsmooth_problem),
          method("topology", &interface<Handle>::topology),
          method("current_step", &interface<Handle>::current_step),
          method("time_step", &interface<Handle>::time_step),
          method("compute_one_step", &interface<Handle>::compute_one_step),
          method("has_next_event", &interface<Handle>::has_next_event),
          method("update_indexsets",
                 &interface<Handle>::update_indexsets<indice>));
    }
  };
};
}  // namespace siconos::simul
