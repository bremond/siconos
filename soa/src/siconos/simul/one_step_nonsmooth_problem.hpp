#pragma once
#include <LinearComplementarityProblem.h>
#include <lcp_cst.h>

#include "SolverOptions.h"
#include "siconos/siconos.hpp"
#include "siconos/utils/pattern.hpp"
#include "siconos/utils/some.hpp"

namespace siconos {
struct lcp : item<> {
  struct options : some::type<pointer<SolverOptions>>, access<options>, some::attribute<> {};
  struct problem : some::type<pointer<LinearComplementarityProblem>>,
                   access<problem>, some::attribute<> {};

  using attributes = gather<options, problem>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) options() { return Handle::type::options::at(*self()); };
    decltype(auto) problem() { return Handle::type::problem::at(*self()); };

    void create(int solver_id)
    {
      self()->options().reset(solver_options_create(solver_id));
      self()->problem().reset(new LinearComplementarityProblem);
    };

    ~interface() { solver_options_delete(&*self()->options()); };
  };
};

template <typename NonsmoothProblem>
struct one_step_nonsmooth_problem : item<> {

  struct problem : some::item_ref<NonsmoothProblem>, access<problem> {};
  struct level : some::indice, access<level> {};

  using attributes = gather<level, problem>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) problem(){return Handle::type::problem::at(*self()); };
    decltype(auto) level() { return Handle ::type ::level ::at(*self()); };

    void solve() { auto& data = self()->data(); }
  };
};

}  // namespace siconos
