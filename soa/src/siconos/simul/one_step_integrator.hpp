#pragma once

#include "siconos/storage/storage.hpp"
#include "siconos/utils/pattern.hpp"

namespace siconos {

template <typename DynamicalSystem, typename Interaction>
struct one_step_integrator {
  using interaction = Interaction;
  using nonsmooth_law = typename interaction::nonsmooth_law;
  using system = DynamicalSystem;
  using y = typename interaction::y;
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
    struct q_vector_assembled : some::unbounded_vector<q>,
                                access<q_vector_assembled> {};
    struct velocity_vector_assembled : some::unbounded_vector<velocity>,
                                       access<velocity_vector_assembled> {};
    struct y_vector_assembled : some::unbounded_vector<y>,
                                access<y_vector_assembled> {};
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
                          y_vector_assembled>;

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
      decltype(auto) y_vector_assembled()
      {
        return Handle ::type ::y_vector_assembled ::at(*self());
      }
      decltype(auto) mass_matrix_assembled()
      {
        return Handle ::type ::mass_matrix_assembled ::at(*self());
      }
      decltype(auto) w_matrix()
      {
        return Handle ::type ::w_matrix ::at(*self());
      }

      bool add_interaction_in_index_set(auto inter, auto h, auto i)
      {
        auto y = inter.y()[i - 1];
        const auto ydot = inter.y()[i];
        const auto &gamma_v = 0.5;
        y += gamma_v * h * ydot;

        return y <= self()->constraint_activation_threshold();
      };

      bool remove_interaction_from_index_set(auto inter, auto h, auto i)
      {
        return !add_interaction_in_index_set(inter, h, i);
      }

      // strategy 1 : assemble the matrix for involved ds only
      auto assemble_h_matrix_for_involved_ds(auto step)
      {
        auto &data = self()->data();
        auto &h_matrices = memory(step, get_memory<h_matrix>(data));
        auto size = std::size(h_matrices);
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
          set_value(self()->h_matrix_assembled(), i, j, mat);
        }
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

      auto assemble_mass_matrix_for_involved_ds(auto step)
      {
        auto &data = self()->data();

        auto size = std::size(memory(step, get_memory<h_matrix>(data)));

        auto &mass_matrices = memory(step, get_memory<mass_matrix>(data));
        auto &involved_ds =
          ground::get<attached_storage<system, symbol<"involved">,
                                       some::boolean>>(data)[0];

        resize(self()->mass_matrix_assembled(), size, size);

        for (auto [i, mat] : ranges::views::enumerate(mass_matrices) |
                                 ranges::views::filter([&involved_ds](auto k_m) {
                                   auto [k, _] = k_m;
                                   return involved_ds[k];
                                 })) {
          set_value(self()->mass_matrix_assembled(), i, i, mat);
          std::cout << size0(mass_matrix_assembled()) << "," << size1(mass_matrix_assembled()) << std::endl;
        }
        assert( size0(mass_matrix_assembled()) == size1(mass_matrix_assembled()));
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

      auto compute_w_matrix(auto step)
      {
        auto &data = self()->data();
        using info_t = std::decay_t<decltype(ground::get<info>(data))>;
        using env = typename info_t::env;
        auto tmp_matrix = typename traits::config<env>::template convert<
            some::unbounded_matrix<some::transposed_matrix<h_matrix>>>::
            type{};

        resize(tmp_matrix, size1(h_matrix_assembled()), size0(h_matrix_assembled()));
        // M^-1 H^t
//        if constexpr (has_property<mass_matrix, some::diagonal>) {
//          prod(inv(mass_matrix_assembled()), trans(h_matrix_assembled(), tmp_matrix));
//        }
//        else  // general case
//        {
        solvet(mass_matrix_assembled(), h_matrix_assembled(), tmp_matrix);
//        }

        // aliasing ?
        resize(w_matrix(), size0(h_matrix_assembled()), size1(tmp_matrix));

        // prod(h_matrix_assembled(), tmp_matrix, w_matrix());
      }

      auto solve_nonsmooth_problem(auto step) {}

      auto compute_iteration_matrix(auto step)
      {
        auto &data = self()->data();
        auto &mass_matrices = get_memory<mass_matrix>(data);
        auto &external_forces = get_memory<fext>(data);

        auto &mats = memory(step, mass_matrices);
        auto &fs = memory(step, external_forces);

        if constexpr (has_property<mass_matrix, some::time_invariant>(data)) {
          if constexpr (has_property<fext, some::time_invariant>(data)) {
            if constexpr (has_property<mass_matrix, property::diagonal>(data)) {
              for (auto [mat, f] : ranges::views::zip(mats, fs)) {
                solve_in_place(mat, f);
              }
            }
          }
        }
      };
      auto compute_free_state(auto step, auto h)
      {
        auto &data = self()->data();
        auto &velocities = get_memory<velocity>(data);
        auto fexts = get_memory<fext>(data);
        auto theta_ = self()->theta();

        auto &vs = memory(step, velocities);
        auto &vs_next = memory(step + 1, velocities);
        auto &minv_fs = memory(step, fexts);
        auto &minv_fs_next = memory(step + 1, fexts);

        // f <- M  f
        for (auto [v, v_next, minv_f, minv_f_next] :
             ranges::views::zip(vs, vs_next, minv_fs, minv_fs_next)) {
          // beware of temporaries...
          v_next = v + h * theta_ * minv_f + h * (1 - theta_) * minv_f_next;
        }
      };
    };
  };
};
}  // namespace siconos
