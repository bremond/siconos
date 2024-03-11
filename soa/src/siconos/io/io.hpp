#pragma once

#include "siconos/storage/pattern/base.hpp"
#include "siconos/storage/pattern/pattern.hpp"
#include "siconos/storage/some/some.hpp"
#include "siconos/storage/storage.hpp"

namespace siconos::io {
using namespace storage;
using namespace storage::pattern;

template <typename System>
struct io : item<> {
  using system = System;
  using attributes =
      gather<attribute<"pos", some::unbounded_collection<some::vector<
                                  some::scalar, some::indice_value<4>>>>,
             attribute<"vel", some::unbounded_collection<some::vector<
                                  some::scalar, some::indice_value<4>>>>>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) positions(auto step)
    {
      auto& data = self()->data();
      using env_t = decltype(self()->env());
      //      using indice = typename env_t::indice;
      using scalar = typename env_t::scalar;

      using system = System;

      auto& ids = storage::prop_values<system, "id">(data, step);
      auto& qs = storage::attr_values<system, "q">(data, step);

      attr<"pos">(*self()).clear();

      for (auto [id, q] : view::zip(ids, qs)) {
        attr<"pos">(*self()).push_back({id, q[0], q[1], q[2]});
      }

      return algebra::matrix_view<algebra::unbounded_matrix<scalar>>(
          attr<"pos">(*self()).data()->data(), attr<"pos">(*self()).size(),
          attr<"pos">(*self()).data()->size() + 1);
    }

    decltype(auto) velocities(auto step)
    {
      auto& data = self()->data();
      using env_t = decltype(self()->env());
      //      using indice = typename env_t::indice;
      using scalar = typename env_t::scalar;

      using system = System;

      auto& ids = storage::prop_values<system, "id">(data, step);
      auto& velos = storage::attr_values<system, "velocities">(data, step);

      attr<"vel">(*self()).clear();

      for (auto [id, velo] : view::zip(ids, velos)) {
        attr<"vel">(*self()).push_back({id, velo[0], velo[1], velo[2]});
      }

      return algebra::matrix_view<algebra::unbounded_matrix<scalar>>(
          attr<"vel">(*self()).data()->data(), attr<"vel">(*self()).size(),
          attr<"vel">(*self()).data()->size() + 1);
    }

    auto methods()
    {
      using env_t = decltype(self()->env());
      using indice = typename env_t::indice;

      return collect(
        method("positions", &interface<Handle>::positions<indice>),
        method("velocities", &interface<Handle>::positions<indice>));
    }
  };
};
}  // namespace siconos::io
