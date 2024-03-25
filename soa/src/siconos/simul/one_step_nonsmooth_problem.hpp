#pragma once

#include <FrictionContactProblem.h>
#include <LinearComplementarityProblem.h>
#include <NonSmoothDrivers.h>
#include <fmt/core.h>
#include <fmt/ranges.h>
#include <lcp_cst.h>
#include <siconos/numerics/Friction_cst.h>

#include "SolverOptions.h"
#include "siconos/simul/simul_head.hpp"

namespace siconos {

namespace simul {

struct solver_options : storage::data_holder<SolverOptions> {
  using attributes = typename storage::data_holder<SolverOptions>::attributes;

  template <typename Handle>
  struct interface : storage::data_holder<SolverOptions>::template interface<
                         Handle> {
    using default_interface<Handle>::self;
    void create(int solver_id = SICONOS_FRICTION_2D_LEMKE)
    {
      self()->instance().reset(solver_options_create(solver_id),
                               [](SolverOptions* so) {
                                 solver_options_delete(so);
                                 delete so;
                               });
    }

    auto methods()
    {
      return collect(method("create", &interface<Handle>::create));
    }
  };
};

template <typename Formulation>
struct nonsmooth_problem : storage::data_holder<Formulation> {
  using attributes = typename storage::data_holder<Formulation>::attributes;

  using formulation_t = Formulation;
  template <typename Handle>
  struct interface : storage::data_holder<Formulation>::template interface<
                         Handle> {
    using default_interface<Handle>::self;

    void create(int solver_id = SICONOS_FRICTION_2D_LEMKE)
    {
      self()->instance().reset(new Formulation);
    };

    //    ~interface() { solver_options_delete(&*_options); };
    auto methods()
    {
      return collect(method("create", &interface<Handle>::create));
    }
  };
};

template <typename NonsmoothProblem>
struct one_step_nonsmooth_problem : item<> {
  using problem_t = NonsmoothProblem;
  using attributes =
      gather<attribute<"level", some::indice>,
             attribute<"options", some::item_ref<solver_options>>,
             attribute<"problem", some::item_ref<NonsmoothProblem>>>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) options()
    {
      return storage::handle(self()->data(), attr<"options">(*self()));
    };
    decltype(auto) problem()
    {
      return storage::handle(self()->data(), attr<"problem">(*self()));
    };
    decltype(auto) level() { return attr<"level">(*self()); };

    template <typename Formulation, match::matrix W, match::vector V>
    void solve(algebra::mat<W>& w_mat, algebra::vec<V>& q_vec,
               algebra::vec<V>& z_vec, algebra::vec<V>& w_vec)
    {
      using fmt::print;
      if constexpr (std::derived_from<Formulation,
                                      LinearComplementarityProblem>) {
        // w_mat cannot be sparse
        auto w_mat_dense = NM_create(NM_DENSE, size0(w_mat), size1(w_mat));
        NM_to_dense(w_mat._m, w_mat_dense);

        /*        print("w_mat_dense:\n");
                  NM_display(w_mat_dense);*/
        self()->problem().instance()->size = size0(w_mat);
        self()->problem().instance()->M = w_mat_dense;
        self()->problem().instance()->q = q_vec._v->matrix0;

        print("LCP [\n");

        print("W:\n");
        algebra::display(w_mat);
        print("----\n");
        print("----\n");

        print("q:\n");
        algebra::display(q_vec);
        print("----\n");

        print("z:\n");
        algebra::display(z_vec);
        print("----\n");

        print("w:\n");
        algebra::display(w_vec);
        print("----\n");

        linearComplementarity_driver(&*self()->problem().instance(),
                                     z_vec._v->matrix0, w_vec._v->matrix0,
                                     &*options().instance());

        print("q:\n");
        algebra::display(q_vec);
        print("----\n");

        print("z:\n");
        algebra::display(z_vec);
        print("----\n");

        print("w:\n");
        algebra::display(w_vec);
        print("----\n");

        print("]\n\n");
        NM_free(w_mat_dense);
      }
      else if constexpr (std::derived_from<Formulation,
                                           FrictionContactProblem>) {
        self()->problem().instance()->dimension = 2;
        self()->problem().instance()->numberOfContacts = size0(w_mat);
        self()->problem().instance()->M = w_mat._m;
        self()->problem().instance()->q = q_vec._v->matrix0;
        self()->problem().instance()->mu =
            (double*)malloc(size0(w_mat) * sizeof(double));

        for (unsigned int i = 0; i < size0(w_mat); ++i)
          self()->problem().instance()->mu[i] = 0.;
        fc2d_driver(&*self()->problem().instance(), z_vec._v->matrix0,
                    w_vec._v->matrix0, &*options().instance());

        free(self()->problem().instance()->mu);
      }
    }

    auto methods()
    {
      return collect(method("options", &interface<Handle>::options),
                     method("problem", &interface<Handle>::problem),
                     method("level", &interface<Handle>::level));
    }
  };
};
}  // namespace simul
}  // namespace siconos
