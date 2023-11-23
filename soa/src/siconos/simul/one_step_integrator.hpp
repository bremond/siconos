#pragma once

#include "siconos/algebra/numerics.hpp"
#include "siconos/model/lagrangian_r.hpp"
#include "siconos/model/model.hpp"
#include "siconos/simul/simul.hpp"
#include "siconos/storage/storage.hpp"
#include "siconos/utils/print.hpp"
#include "siconos/utils/range.hpp"
#include "siconos/utils/variant.hpp"

namespace siconos::simul {

template <typename DynamicalSystem, typename Interaction>
struct one_step_integrator {
  using interaction = Interaction;
  using nonsmooth_law = typename interaction::nslaw;
  using system = DynamicalSystem;
  using dof = typename interaction::dof;
  using nslaw_size = typename interaction::nslaw_size;
  using nslaw = typename interaction::nslaw;
  using y = attr_t<interaction, "y">;
  using ydot = attr_t<interaction, "ydot">;
  using lambda = attr_t<interaction, "lambda">;
  using relation = attr_t<interaction, "relation">;
  using h_matrix1 = attr_t<interaction, "h_matrix1">;
  using h_matrix2 = attr_t<interaction, "h_matrix2">;

  using q = attr_t<system, "q">;
  using velocity = attr_t<system, "velocity">;
  using fext = attr_t<system, "fext">;

  struct euler : item<> {
    using properties = gather<storage::keep<attr_t<system, "q">, 2>,
                              storage::keep<attr_t<system, "velocity">, 2>>;

    using attributes = gather<>;

    template <typename Handle>
    struct interface : default_interface<Handle> {
      using default_interface<Handle>::self;
      void compute_free_state(auto step, auto h){
          // auto& data = self()->data();
          // auto& velocities = storage::get_memory<velocity>(data);
          // auto& mass_matrices = storage::get_memory<mass_matrix>(data);
          // auto& external_forces = storage::get_memory<fext>(data);

          // auto& Ms =      storage::memory(step, mass_matrices);
          // auto& vs =      storage::memory(step, velocities);
          // auto& vs_next = storage::memory(step+1, velocities);
          // auto& fs =      storage::memory(step, external_forces);

      };
    };
  };

  struct moreau_jean : item<> {

    using attributes = types::attributes<
        attribute<"theta", some::scalar>, attribute<"gamma", some::scalar>,
        attribute<"constraint_activation_threshold", some::scalar>,
        attribute<"h_matrix_assembled", some::unbounded_matrix<h_matrix1>>,
        attribute<"mass_matrix_assembled",
                  some::unbounded_matrix<attr_t<system, "mass_matrix">>>,
        attribute<"w_matrix",
                   some::unbounded_matrix<some::matrix<
                       some::scalar, nth_t<0, typename h_matrix1::sizes>,
                       nth_t<0, typename h_matrix1::sizes>>>>,
        attribute<
            "q_nsp_vector_assembled",
            some::unbounded_vector<some::vector<some::scalar, nslaw_size>>>,
        attribute<"velocity_vector_assembled",
                  some::unbounded_vector<velocity>>,
        attribute<"y_vector_assembled", some::unbounded_vector<ydot>>,
        attribute<"ydot_vector_assembled", some::unbounded_vector<ydot>>,
        attribute<"free_velocity_vector_assembled",
                  some::unbounded_vector<attr_t<system, "velocity">>>,
        attribute<"lambda_vector_assembled", some::unbounded_vector<lambda>>,
        attribute<"p0_vector_assembled",
                  some::unbounded_vector<some::vector<some::scalar, dof>>>>;

    using properties = gather<storage::keep<attr_t<system, "q">, 2>,
                              storage::keep<attr_t<system, "velocity">, 2>,
                              storage::keep<y, 2>, storage::keep<ydot, 2>>;

    template <typename Handle>
    struct interface : default_interface<Handle> {
      using default_interface<Handle>::self;

      decltype(auto) theta() { return attr<"theta">(*self()); }
      decltype(auto) gamma() { return attr<"gamma">(*self()); }
      decltype(auto) constraint_activation_threshold()
      {
        return attr<"constraint_activation_threshold">(*self());
      }
      decltype(auto) h_matrix_assembled()
      {
        return attr<"h_matrix_assembled">(*self());
      }
      decltype(auto) q_nsp_vector_assembled()
      {
        return attr<"q_nsp_vector_assembled">(*self());
      }
      decltype(auto) velocity_vector_assembled()
      {
        return attr<"velocity_vector_assembled">(*self());
      }
      decltype(auto) free_velocity_vector_assembled()
      {
        return attr<"free_velocity_vector_assembled">(*self());
      }
      decltype(auto) lambda_vector_assembled()
      {
        return attr<"lambda_vector_assembled">(*self());
      }
      decltype(auto) p0_vector_assembled()
      {
        return attr<"p0_vector_assembled">(*self());
      }
      decltype(auto) y_vector_assembled()
      {
        return attr<"y_vector_assembled">(*self());
      }
      decltype(auto) ydot_vector_assembled()
      {
        return attr<"ydot_vector_assembled">(*self());
      }
      decltype(auto) mass_matrix_assembled()
      {
        return attr<"mass_matrix_assembled">(*self());
      }
      decltype(auto) w_matrix()
      {
        return attr<"w_matrix">(*self());
      }

      void compute_output(auto step)
      {
        auto &data = self()->data();

        auto &ys = storage::attr_values<y>(data, step);
        auto &ydots = storage::attr_values<ydot>(data, step);
        auto &h_matrices1 = storage::attr_values<h_matrix1>(data, step);
        auto &h_matrices2 = storage::attr_values<h_matrix2>(data, step);

        auto &qs = storage::attr_values<system, "q">(data, step);
        auto &velocities =
            storage::attr_values<attr_t<system, "velocity">>(data, step);

        auto &ds1s = storage::prop_values<interaction, "ds1">(data, step);
        auto &ds2s = storage::prop_values<interaction, "ds2">(data, step);

        auto &ndss = storage::prop_values<interaction, "nds">(data, step);

        // global h_matrix is not assembled at this stage
        for (auto [y, ydot, hm1, hm2, ds1, ds2, nds] : views::zip(
                 ys, ydots, h_matrices1, h_matrices2, ds1s, ds2s, ndss)) {
          y = hm1 * qs[ds1.get()];
          ydot = hm1 * velocities[ds1.get()];

          if (nds == 2) {
            y += hm2 * qs[ds2.get()];
            ydot += hm2 * velocities[ds2.get()];
          }
        }
      }

      void compute_input()
      {
        auto &h_matrix = h_matrix_assembled();
        auto &lambda = lambda_vector_assembled();
        auto &ydot = ydot_vector_assembled();
        auto &p0 = p0_vector_assembled();
        auto &velo = velocity_vector_assembled();
        auto &mass_matrix = mass_matrix_assembled();

        resize(p0, size1(h_matrix));
        resize(velo, size1(h_matrix));

        transpose(h_matrix);

        // p0 += h_matrix^t * lambda
        prodt1(h_matrix, lambda, p0);

        // velo <- mass_matrix^-1 * p0
        solve(mass_matrix, p0, velo);

        // velo += h_matrix^t * ydot
        prodt1(h_matrix, ydot, velo);

        /*      print("ydot assembled:\n");
                numerics::display(ydot);
                print("lambda assembled:\n");
                numerics::display(lambda);
                print("p0 assembled:\n");
                numerics::display(p0);
                print("velo assembled:\n");
                numerics::display(velo);*/
      }

      auto compute_active_interactions(auto step, auto h)
      {
        auto &data = self()->data();

        using info_t =
            std::decay_t<decltype(ground::get<storage::info>(data))>;
        using env = typename info_t::env;
        using indice = typename env::indice;

        auto &ys = storage::attr_values<y>(data, step + 1);
        auto &ydots = storage::attr_values<ydot>(data, step + 1);

        auto &ids1s = storage::prop_values<interaction, "ds1">(data, step);
        auto &ids2s = storage::prop_values<interaction, "ds2">(data, step);

        auto &activations =
            storage::prop_values<interaction, "activation">(data, step);

        auto &involveds =
            storage::prop_values<system, "involved">(data, step);

        const auto &interactions = storage::handles<interaction>(data, step);

        // all ds -> not involved
        // without zip : involved is a copy not a ref!!
        for (auto [involved] : views::zip(involveds)) {
          involved = false;
        };

        auto gamma_v = 0.5;

        indice ds_counter = 0;
        indice inter_counter = 0;
        for (auto [y, ydot, activation, ids1, ids2, inter] :
             views::zip(ys, ydots, activations, ids1s, ids2s, interactions)) {
          auto b = siconos::variant::visit(
              data, inter.relation(),
              ground::overload(
                  [](match::linear_relation auto &real_relation) {
                    return real_relation.b();
                  },
                  // no b() present
                  [](auto) { return 0.; }));
          // on normal component
          activation = ((y + gamma_v * h * ydot)(0) + b <=
                        self()->constraint_activation_threshold());

          if (activation) {
            inter_counter++;

            auto ds1 = storage::handle(data, ids1);
            auto ds2 = storage::handle(data, ids2);

            if (!prop<"involved">(ds1)) {
              prop<"involved">(ds1) = true;
              prop<"index">(ds1) = ds_counter++;
            };

            if (!prop<"involved">(ds2)) {
              prop<"involved">(ds2) = true;
              prop<"index">(ds2) = ds_counter++;
            }

            //          print(
            //              "\nstep: {}, time: {} => ACTIVATION {}<->{}
            //              !\ny:{}, " "ydot:{}\n",
            //             step, step * h, ids1.get(), ids2.get(), y, ydot);
          }
        }
        return std::pair{inter_counter, ds_counter};
      }

      // strategy 1 : assemble the matrix for involved ds only
      auto assemble_h_matrix_for_involved_ds(auto step, auto ninter, auto nds)
      {
        auto &data = self()->data();

        resize(self()->h_matrix_assembled(), ninter, nds);

        for (auto [i, hi] :
             (storage::handles<interaction>(data, step) |
              views::filter([](auto h) { return prop<"activation">(h); })) |
                 views::enumerate) {
          auto mat1 = hi.h_matrix1();
          auto mat2 = hi.h_matrix2();
          auto ids1 = prop<"ds1">(hi);
          auto ids2 = prop<"ds2">(hi);

          auto j1 = prop<"index">(storage::handle(data, ids1));
          auto j2 = prop<"index">(storage::handle(data, ids2));

          if (j1 == j2) {
            // one block
            set_value(self()->h_matrix_assembled(), i, j1,
                      mat1);  // sparse block matrix
          }
          else {
            // i!=j blocks
            set_value(self()->h_matrix_assembled(), i, j1,
                      mat1);  // sparse block matrix
            set_value(self()->h_matrix_assembled(), i, j2,
                      mat2);  // sparse block matrix
          }
        }

        /*      print("h_matrix:\n");
                numerics::display(h_matrix_assembled());
                print("================\n");*/
      }

      auto resize_assembled_vectors(auto step, auto ninter)
      {
        resize(self()->y_vector_assembled(), ninter);
        resize(self()->ydot_vector_assembled(), ninter);
        resize(self()->lambda_vector_assembled(), ninter);
      }

      auto assemble_mass_matrix_for_involved_ds(auto step, auto size)
      {
        auto &data = self()->data();

        // size may be 0
        resize(mass_matrix_assembled(), size, size);
        resize(free_velocity_vector_assembled(), size);  // !!!

        for (auto hds :
             storage::handles<system>(data, step) |
                 views::filter([](auto h) { return prop<"involved">(h); })) {
          auto i = prop<"index">(hds);
          set_value(mass_matrix_assembled(), i, i, hds.mass_matrix());
        }

        /*      print("mass_matrix:\n");
                numerics::display(mass_matrix_assembled());
                print("================\n");
                assert(size0(mass_matrix_assembled()) ==
                size1(mass_matrix_assembled()));*/
      }

      // compute H vfree
      auto compute_q_nsp_vector_assembled(auto step, auto ninter)
      {
        auto &data = self()->data();

        resize(q_nsp_vector_assembled(), ninter);

        auto &ydots = storage::attr_values<ydot>(data, step + 1);

        auto k = 0;
        for (auto [i, inter] : storage::handles<interaction>(data, step + 1) |
                                   views::enumerate) {
          if (prop<"activation">(inter)) {
            set_value(q_nsp_vector_assembled(), k++, ydots[inter.get()]);
          };
        }
      }
      // compute H M^-1 H^t
      auto compute_w_matrix(auto step)
      {
        auto &data = self()->data();
        using info_t =
            std::decay_t<decltype(ground::get<storage::info>(data))>;
        using env = typename info_t::env;
        auto tmp_matrix = typename traits::config<env>::template convert<
            some::unbounded_matrix<some::transposed_matrix<h_matrix1>>>::
            type{};

        resize(tmp_matrix, size1(h_matrix_assembled()),
               size0(h_matrix_assembled()));
        // M^-1 H^t
        //        if constexpr (has_property<mass_matrix, some::diagonal>) {
        //          prod(inv(mass_matrix_assembled()),
        //          trans(h_matrix_assembled(), tmp_matrix));
        //        }
        //        else  // general case
        //        {
        solvet(mass_matrix_assembled(), h_matrix_assembled(), tmp_matrix);
        //        }

        // aliasing ?
        resize(w_matrix(), size0(h_matrix_assembled()), size1(tmp_matrix));

        prod(h_matrix_assembled(), tmp_matrix, w_matrix());

        /*      print("w_matrix:\n");
                numerics::display(w_matrix());
                print("================\n");*/
      }

      void update_velocity_for_involved_ds() {}

      void nsl_effect_on_free_output(auto step)
      {
        auto &data = self()->data();
        auto &ydots = storage::attr_values<ydot>(data, step);
        auto &ydots_next = storage::attr_values<ydot>(data, step + 1);

        auto &inslaws =
            storage::attr_values<attr_t<interaction, "nslaw">>(data, step);

        for (auto [ydot, ydot_next, inslaw] :
             views::zip(ydots, ydots_next, inslaws)) {
          ydot_next += storage::handle(data, inslaw).e() * ydot;
        }
      }

      void update_all_velocities(auto step)
      {
        auto &data = self()->data();
        auto &velo = velocity_vector_assembled();
        auto &all_vs = storage::attr_values<velocity>(data, step + 1);
        auto &involved_ds =
            storage::prop_values<system, "involved">(data, step);

        auto &indices = storage::prop_values<system, "index">(data, step);

        // involved ds velocities -> ds velocities
        for (auto [i, iv] : views::zip(indices, all_vs) | views::enumerate |
                                views::filter([&involved_ds](auto k__) {
                                  auto [k, _] = k__;
                                  return involved_ds[k];
                                })) {
          auto [indx, v] = iv;
          // copy
          v += get_vector(velo, indx);
        }
      }

      auto update_positions(auto step, auto h)
      {
        auto &data = self()->data();

        auto &xs = storage::attr_values<system, "q">(data, step);
        auto &xs_next = storage::attr_values<system, "q">(data, step + 1);
        auto &vs = storage::attr_values<velocity>(data, step);
        auto &vs_next = storage::attr_values<velocity>(data, step + 1);

        for (auto [x, x_next, v, v_next] :
             views::zip(xs, xs_next, vs, vs_next)) {
          x_next = x + h * (theta() * v + (1.0 - theta()) * v_next);
        }
      }

      auto compute_iteration_matrix(auto step)
      {
        auto &data = self()->data();
        auto &mass_matrices =
            storage::attr_memory<system, "mass_matrix">(data);
        auto &external_forces = storage::attr_memory<fext>(data);

        auto &mats = storage::memory(step, mass_matrices);
        auto &fs = storage::memory(step, external_forces);

        for (auto [mat, f] : views::zip(mats, fs)) {
          f = mat.inverse() * f;
        }
      }

      auto compute_free_state(auto step, auto h)
      {
        auto &data = self()->data();
        auto &velocities = storage::attr_memory<velocity>(data);
        auto &fexts = storage::attr_memory<fext>(data);
        auto &theta_ = self()->theta();

        auto &vs = storage::memory(step, velocities);
        auto &vs_next = storage::memory(step + 1, velocities);
        auto &minv_fs = storage::memory(step, fexts);
        auto &minv_fs_next = storage::memory(step + 1, fexts);

        // for all ds
        for (auto [v, v_next, minv_f, minv_f_next] :
             views::zip(vs, vs_next, minv_fs, minv_fs_next)) {
          // note: theta useless if fext is constant
          v_next = v + h * theta_ * minv_f + h * (1 - theta_) * minv_f_next;
        }
      }

      void initialize(auto step)
      {
        compute_iteration_matrix(step);
        compute_h_matrices(step);
      }

      void compute_h_matrices(auto step)
      {
        auto &data = self()->data();

        auto &h_matrices1 = storage::attr_values<h_matrix1>(data, step);
        auto &h_matrices2 = storage::attr_values<h_matrix2>(data, step);
        auto &ndss = storage::prop_values<interaction, "nds">(data, step);

        auto &ds1s = storage::prop_values<interaction, "ds1">(data, step);
        auto &ds2s = storage::prop_values<interaction, "ds2">(data, step);
        auto &relations = storage::attr_values<relation>(data, step);

        for (auto [rel, hm1, hm2, nds, ds1, ds2] : views::zip(
                 relations, h_matrices1, h_matrices2, ndss, ds1s, ds2s)) {
          // local binding not enough to be passed to lambda...
          auto &hhm1 = hm1;
          auto &hhm2 = hm2;
          auto &q1 = storage::handle(data, ds1).q();

          if (nds == 2) {
            auto &q2 = storage::handle(data, ds2).q();

            siconos::variant::visit(data, rel,
                                    ground::overload(
                                        [&step, &hhm1, &hhm2, &q1,
                                         &q2](match::relation2 auto &rrel) {
                                          rrel.compute_jachq(step, q1, q2,
                                                             hhm1, hhm2);
                                        },
                                        [](auto) { assert(false); }));
          }
          else {
            // nds == 1
            siconos::variant::visit(
                data, rel,
                ground::overload(
                    [&step, &hhm1, &q1](match::relation1 auto rrel) {
                      rrel.compute_jachq(step, q1, hhm1);
                    },
                    [](auto rrel) { assert(false); }));
          }
        };
      }

      void update_h_matrices(auto step)
      {
        auto &data = self()->data();
        using data_t = const std::decay_t<decltype(data)>;
        if constexpr (!storage::has_property_t<
                          interaction, property::time_invariant, data_t>()) {
          compute_h_matrices(step);
        }
      };

      void update_iteration_matrix(auto current_step)
      {
        using data_t = const std::decay_t<decltype(self()->data())>;
        if constexpr (!storage::has_property_t<typename system::fext,
                                               property::time_invariant,
                                               data_t>()) {
          // constant fext => constant iteration matrix
          compute_iteration_matrix(current_step);
        }
      }
    };
  };
};
}  // namespace siconos::simul
