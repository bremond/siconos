#pragma once

#include <fmt/core.h>
#include <fmt/ranges.h>

#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/filter.hpp>
#include <range/v3/view/zip.hpp>

#include "siconos/storage/storage.hpp"
#include "siconos/utils/pattern.hpp"
using fmt::print;

namespace siconos {

template <typename DynamicalSystem, typename Interaction>
struct one_step_integrator {
  using interaction = Interaction;
  using nonsmooth_law = typename interaction::nonsmooth_law;
  using system = DynamicalSystem;
  using dof = typename interaction::dof;
  using nslaw_size = typename interaction::nslaw_size;
  using y = typename interaction::y;
  using ydot = typename interaction::ydot;
  using lambda = typename interaction::lambda;
  using q = typename system::q;
  using velocity = typename system::velocity;
  using mass_matrix = typename system::mass_matrix;
  using h_matrix = typename interaction::h_matrix;
  using fext = typename system::fext;

  struct euler : item<> {
    using properties = gather<keep<q, 2>, keep<velocity, 2>>;

    using attributes = gather<>;

    template <typename Handle>
    struct interface : default_interface<Handle> {
      using default_interface<Handle>::self;
      void compute_free_state(auto step, auto h){
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

  struct moreau_jean : item<> {
    struct theta : some::scalar, access<theta> {};
    struct gamma : some::scalar, access<gamma> {};
    struct constraint_activation_threshold
        : some::scalar,
          access<constraint_activation_threshold> {};

    struct h_matrix_assembled : some::unbounded_matrix<h_matrix>,
                                access<h_matrix_assembled> {};
    struct q_vector_assembled
        : some::unbounded_vector<some::vector<some::scalar, nslaw_size>>,
          access<q_vector_assembled> {};
    struct velocity_vector_assembled : some::unbounded_vector<velocity>,
                                       access<velocity_vector_assembled> {};
    struct free_velocity_vector_assembled
        : some::unbounded_vector<velocity>,
          access<free_velocity_vector_assembled> {};
    struct y_vector_assembled : some::unbounded_vector<y>,
                                access<y_vector_assembled> {};

    struct ydot_vector_assembled : some::unbounded_vector<ydot>,
                                   access<ydot_vector_assembled> {};

    struct lambda_vector_assembled : some::unbounded_vector<lambda>,
                                     access<lambda_vector_assembled> {};

    struct p0_vector_assembled
        : some::unbounded_vector<some::vector<some::scalar, dof>>,
          access<p0_vector_assembled> {};
    struct mass_matrix_assembled : some::unbounded_matrix<mass_matrix>,
                                   access<mass_matrix_assembled> {};
    struct w_matrix
        : some::unbounded_matrix<
              some::matrix<some::scalar, nth_t<0, typename h_matrix::sizes>,
                           nth_t<0, typename h_matrix::sizes>>>,
          access<w_matrix> {};

    using attributes =
        types::attributes<theta, gamma, constraint_activation_threshold,
                          h_matrix_assembled, mass_matrix_assembled, w_matrix,
                          q_vector_assembled, velocity_vector_assembled,
                          y_vector_assembled, ydot_vector_assembled,
                          free_velocity_vector_assembled,
                          lambda_vector_assembled, p0_vector_assembled>;

    using properties = gather<keep<typename system::q, 2>,
                              keep<typename system::velocity, 2>>;

    template <typename Handle>
    struct interface : default_interface<Handle> {
      using default_interface<Handle>::self;

      decltype(auto) theta() { return Handle ::type ::theta ::at(*self()); }
      decltype(auto) gamma() { return Handle ::type ::gamma ::at(*self()); }
      decltype(auto) constraint_activation_threshold()
      {
        return Handle ::type ::constraint_activation_threshold ::at(*self());
      }
      decltype(auto) h_matrix_assembled()
      {
        return Handle ::type ::h_matrix_assembled ::at(*self());
      }
      decltype(auto) q_vector_assembled()
      {
        return Handle ::type ::q_vector_assembled ::at(*self());
      }
      decltype(auto) velocity_vector_assembled()
      {
        return Handle ::type ::velocity_vector_assembled ::at(*self());
      }
      decltype(auto) free_velocity_vector_assembled()
      {
        return Handle ::type ::free_velocity_vector_assembled ::at(*self());
      }
      decltype(auto) lambda_vector_assembled()
      {
        return Handle ::type ::lambda_vector_assembled ::at(*self());
      }
      decltype(auto) p0_vector_assembled()
      {
        return Handle ::type ::p0_vector_assembled ::at(*self());
      }
      decltype(auto) y_vector_assembled()
      {
        return Handle ::type ::y_vector_assembled ::at(*self());
      }
      decltype(auto) ydot_vector_assembled()
      {
        return Handle ::type ::ydot_vector_assembled ::at(*self());
      }
      decltype(auto) mass_matrix_assembled()
      {
        return Handle ::type ::mass_matrix_assembled ::at(*self());
      }
      decltype(auto) w_matrix()
      {
        return Handle ::type ::w_matrix ::at(*self());
      }

      void compute_output(auto step)
      {
        auto &data = self()->data();

        auto &ys = memory(step, get_memory<y>(data));
        auto &ydots = memory(step, get_memory<ydot>(data));
        auto &h_matrices = memory(step, get_memory<h_matrix>(data));

        auto &qs = memory(step, get_memory<q>(data));
        auto &velocities = memory(step, get_memory<velocity>(data));

        auto &ids1s =
            ground::get<attached_storage<interaction, symbol<"ds1">,
                                         some::item_ref<system>>>(data)[0];
        auto &ids2s =
            ground::get<attached_storage<interaction, symbol<"ds2">,
                                         some::item_ref<system>>>(data)[0];

        for (auto [y, ydot, hm, ids1, jds2] :
             ranges::views::zip(ys, ydots, h_matrices, ids1s, ids2s)) {
          auto i = handle(ids1, data).property(symbol<"index">{});
          y = hm * qs[i];
          ydot = hm * velocities[i];
          // fix for second ds
        }
      }
      void compute_active_interactions(auto step, auto h)
      {
        auto &data = self()->data();

        auto &ys = memory(step, get_memory<y>(data));
        auto &ydots = memory(step, get_memory<ydot>(data));

        auto &activations =
            ground::get<attached_storage<interaction, symbol<"activation">,
                                         some::boolean>>(data)[0];

        auto gamma_v = 0.5;

        for (auto [y, ydot, activation] :
             ranges::views::zip(ys, ydots, activations)) {
          // on normal component
          activation = ((y + gamma_v * h * ydot)(0) <=
                        self()->constraint_activation_threshold());
        }
      }

      // strategy 1 : assemble the matrix for involved ds only
      auto assemble_h_matrix_for_involved_ds(auto step, auto size)
      {
        auto &data = self()->data();
        auto &h_matrices = memory(step, get_memory<h_matrix>(data));
        auto &ids1s =
            ground::get<attached_storage<interaction, symbol<"ds1">,
                                         some::item_ref<system>>>(data)[0];
        auto &ids2s =
            ground::get<attached_storage<interaction, symbol<"ds2">,
                                         some::item_ref<system>>>(data)[0];

        resize(self()->h_matrix_assembled(), size, size);
        for (auto [mat, ids1, ids2] :
             ranges::views::zip(h_matrices, ids1s, ids2s)) {
          auto i = handle(ids1, data)
                       .property(symbol<"index">{});  // involved index
          auto j = handle(ids2, data)
                       .property(symbol<"index">{});  // involved index
          set_value(self()->h_matrix_assembled(), i, j,
                    mat);  // sparse block matrix
        }
      }

      auto resize_assembled_vectors(auto step)
      {
        auto &data = self()->data();
        auto &h_matrices = memory(step, get_memory<h_matrix>(data));
        auto size = std::size(h_matrices);

        resize(self()->y_vector_assembled(), size);
        resize(self()->ydot_vector_assembled(), size);
        resize(self()->lambda_vector_assembled(), size);
      }

      // strategy 2 : assemble the whole matrix (size = number of ds)
      auto assemble_h_matrix_for_all_ds(auto step)
      {
        auto &data = self()->data();

        // size is the number of ds
        auto size = std::size(memory(step, get_memory<mass_matrix>(data)));

        auto &h_matrices = memory(step, get_memory<h_matrix>(data));
        auto &ids1s =
            ground::get<attached_storage<interaction, symbol<"ds1">,
                                         some::item_ref<system>>>(data)[0];
        auto &ids2s =
            ground::get<attached_storage<interaction, symbol<"ds2">,
                                         some::item_ref<system>>>(data)[0];

        resize(self()->h_matrix_assembled(), size, size);
        for (auto [mat, ids1, ids2] :
             ranges::views::zip(h_matrices, ids1s, ids2s)) {
          auto i = handle(ids1, data).get();  // global index
          auto j = handle(ids2, data).get();  // global index
          set_value(self()->h_matrix_assembled(), i, j, mat);
        }
      }

      auto assemble_mass_matrix_for_involved_ds(auto step, auto size)
      {
        auto &data = self()->data();

        auto &mass_matrices = memory(step, get_memory<mass_matrix>(data));
        auto &free_velocities = memory(step + 1, get_memory<velocity>(data));
        auto &involved_ds = ground::get<
            attached_storage<system, symbol<"involved">, some::boolean>>(
            data)[0];

        // size may be 0
        resize(self()->mass_matrix_assembled(), size, size);
        resize(self()->free_velocity_vector_assembled(), size);

        for (auto [i, mat_velo] :
             ranges::views::zip(mass_matrices, free_velocities) |
                 ranges::views::enumerate |
                 ranges::views::filter([&involved_ds](auto k_m) {
                   auto [k, _] = k_m;
                   return involved_ds[k];
                 })) {
          auto &[mat, velo] = mat_velo;
          print("velo {}\n", velo);
          set_value(self()->mass_matrix_assembled(), i, i, mat);
          set_value(self()->free_velocity_vector_assembled(), i, velo);
        }
        assert(size0(mass_matrix_assembled()) ==
               size1(mass_matrix_assembled()));
      }

      auto assemble_mass_matrix_for_all_ds(auto step)
      {
        auto &data = self()->data();

        // size is the number of ds
        auto size = std::size(memory(step, get_memory<mass_matrix>(data)));

        auto &mass_matrices = memory(step, get_memory<mass_matrix>(data));

        self()->mass_matrix_assembled().resize(size, size);

        // also mapping block vector -> block diagonal matrix
        for (auto [i, mat] : ranges::views::enumerate(mass_matrices)) {
          set_value(self()->mass_matrix_assembled(), i, i, mat);
        }
      }

      // compute H vfree
      auto compute_q_vector_assembled(auto step)
      {
        resize(self()->q_vector_assembled(),
               size0(self()->free_velocity_vector_assembled()));
        prod(self()->h_matrix_assembled(),
             self()->free_velocity_vector_assembled(),
             self()->q_vector_assembled());
      }
      // compute H M^-1 H^t
      auto compute_w_matrix(auto step)
      {
        auto &data = self()->data();
        using info_t = std::decay_t<decltype(ground::get<info>(data))>;
        using env = typename info_t::env;
        auto tmp_matrix = typename traits::config<env>::template convert<
            some::unbounded_matrix<some::transposed_matrix<h_matrix>>>::
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
      }

      void update_velocity_for_involved_ds()
      {
        auto &p0 = self()->p0_vector_assembled();
        auto &velo = self()->velocity_vector_assembled();
        auto &mass_matrix = self()->mass_matrix_assembled();

        resize(velo, size0(p0));
        solve(mass_matrix, p0, velo);
      }

      void apply_nonsmooth_law_effect()
      {
        // auto &data = self()->data();
        // auto &velo = self()->velocity_vector_assembled();
        // auto &all_vs = memory(step, get_memory<velocity>(data));
        // auto &all_next_vs = memory(step + 1, get_memory<velocity>(data));
        // auto &nslaws = memory(step, get_memory<typename interaction::nonsmooth_law>(data));
        // auto &involved_ds = ground::get<
        //     attached_storage<system, symbol<"involved">, some::boolean>>(
        //       data)[0];

        // // involved ds velocities -> ds velocities
        // for (auto [i, v_vn] :
        //        ranges::view::zip(all_vs, all_next_vs) | ranges::views::enumerate |
        //        ranges::views::filter([&involved_ds](auto k_m) {
        //            auto [k, _] = k_m;
        //            return involved_ds[k];
        //        })) {
        //   // copy
        //   v += get_vector(velo, i);
        // }

      }

      void update_all_velocities(auto step)
      {
        auto &data = self()->data();
        auto &velo = self()->velocity_vector_assembled();
        auto &all_vs = memory(step + 1, get_memory<velocity>(data));
        auto &involved_ds = ground::get<
            attached_storage<system, symbol<"involved">, some::boolean>>(
            data)[0];

        // involved ds velocities -> ds velocities
        for (auto [i, v] :
             all_vs | ranges::views::enumerate |
                 ranges::views::filter([&involved_ds](auto k_m) {
                   auto [k, _] = k_m;
                   return involved_ds[k];
                 })) {
          // copy
          v += get_vector(velo, i);
        }
      }

      auto update_positions(auto step, auto h)
      {
        auto &data = self()->data();
        auto &xs = memory(step, get_memory<q>(data));
        auto &xs_next = memory(step + 1, get_memory<q>(data));
        auto &vs = memory(step + 1, get_memory<velocity>(data));

        for (auto [x, x_next, v] : ranges::views::zip(xs, xs_next, vs)) {
          x_next = x + h * v;
        }
      }

      auto compute_iteration_matrix(auto step)
      {
        auto &data = self()->data();
        auto &mass_matrices = get_memory<mass_matrix>(data);
        auto &external_forces = get_memory<fext>(data);

        auto &mats = memory(step, mass_matrices);
        auto &fs = memory(step, external_forces);

        //        if constexpr (has_property<mass_matrix,
        //        some::time_invariant>(data)) {
        //          if constexpr (has_property<fext,
        //          some::time_invariant>(data)) {
        //            if constexpr (has_property<mass_matrix,
        //            property::diagonal>(
        //                              data)) {
        for (auto [mat, f] : ranges::views::zip(mats, fs)) {
          f = mat.inverse() * f;
        }
        //            }
        //          }
        //      }
      }

      auto compute_free_state(auto step, auto h)
      {
        auto &data = self()->data();
        auto &velocities = get_memory<velocity>(data);
        auto &fexts = get_memory<fext>(data);
        auto &theta_ = self()->theta();

        auto &vs = memory(step, velocities);
        auto &vs_next = memory(step + 1, velocities);
        auto &minv_fs = memory(step, fexts);
        auto &minv_fs_next = memory(step + 1, fexts);

        // for all ds
        for (auto [v, v_next, minv_f, minv_f_next] :
             ranges::views::zip(vs, vs_next, minv_fs, minv_fs_next)) {
          // theta useless if fext is constant
          v_next = v + h * theta_ * minv_f + h * (1 - theta_) * minv_f_next;
        }
      };
    };
  };
};
}  // namespace siconos
