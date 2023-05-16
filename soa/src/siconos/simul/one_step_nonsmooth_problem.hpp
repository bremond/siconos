#pragma once
#include <LinearComplementarityProblem.h>
#include <NonSmoothDrivers.h>
#include <lcp_cst.h>

#include "SolverOptions.h"
#include "siconos/siconos.hpp"
#include "siconos/utils/pattern.hpp"
#include "siconos/utils/some.hpp"

namespace siconos {

namespace numerics {

template <typename Struct>
struct data_holder : item<> {
  struct instance : some::specific<pointer<Struct>>,
                    access<instance>,
                    some::attribute<> {};

  using attributes = gather<instance>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) instance() { return Handle::type::instance::at(*self()); };
  };
};

struct solver_options : data_holder<SolverOptions> {
  using attributes = typename data_holder<SolverOptions>::attributes;

  template <typename Handle>
  struct interface : data_holder<SolverOptions>::template interface<Handle> {
    using default_interface<Handle>::self;
    void create(int solver_id = SICONOS_LCP_LEMKE)
    {
      self()->instance().reset(solver_options_create(solver_id),
                               [](SolverOptions* so) {
                                 solver_options_delete(so);
                                 delete so;
                               });
    }
  };
};

template <typename Formulation>
struct nonsmooth_problem : data_holder<Formulation> {
  using attributes = typename data_holder<Formulation>::attributes;

  template <typename Handle>
  struct interface : data_holder<Formulation>::template interface<Handle> {
    using default_interface<Handle>::self;

    void create(int solver_id = SICONOS_LCP_LEMKE)
    {
      self()->instance().reset(new Formulation);
    };

    //    ~interface() { solver_options_delete(&*_options); };
  };
};

template <typename NonsmoothProblem>
struct one_step_nonsmooth_problem : item<> {
  struct options : some::item_ref<solver_options>, access<options> {};
  struct problem : some::item_ref<NonsmoothProblem>, access<problem> {};
  struct level : some::indice, access<level> {};

  using attributes = gather<level, options, problem>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) options()
    {
      return handle(Handle::type::options::at(*self()), self()->data());
    };
    decltype(auto) problem()
    {
      return handle(Handle::type::problem::at(*self()), self()->data());
    };
    decltype(auto) level() { return Handle ::type ::level ::at(*self()); };

    template <typename Formulation, siconos::match::matrix W,
              siconos::match::vector V>
    void solve(mat<W>& w_mat, vec<V>& q_vec, vec<V>& z_vec, vec<V>& w_vec)
    {
      if constexpr (std::derived_from<Formulation,
                                      LinearComplementarityProblem>) {
        NM_to_dense(w_mat._m, w_mat._m);
        self()->problem().instance()->size = size0(w_mat);
        self()->problem().instance()->M = w_mat._m;
        self()->problem().instance()->q = q_vec._v->matrix0;

        linearComplementarity_driver(&*self()->problem().instance(),
                                     z_vec._v->matrix0, w_vec._v->matrix0,
                                     &*options().instance());
      }
    }
  };
};
}  // namespace numerics
}  // namespace siconos
