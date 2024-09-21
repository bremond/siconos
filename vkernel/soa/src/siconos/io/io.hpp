#pragma once

#include <stdexcept>

#include "siconos/algebra/numerics.hpp"
#include "siconos/collision/collision.hpp"
#include "siconos/collision/diskfdisk_r.hpp"
#include "siconos/collision/disksegment_r.hpp"
#include "siconos/storage/pattern/base.hpp"
#include "siconos/storage/pattern/pattern.hpp"
#include "siconos/storage/some/some.hpp"
#include "siconos/storage/storage.hpp"
#include "siconos/utils/variant.hpp"

namespace siconos::io {
using namespace storage;
using namespace storage::pattern;

template <typename Osi>
struct io : item<> {
  using osi = Osi;
  using system = typename osi::system;
  using interaction = typename osi::interaction;

  using attributes = gather<
      attribute<"pos_info", some::unbounded_collection<some::vector<
                                some::scalar, some::indice_value<4>>>>,
      attribute<"vel_info", some::unbounded_collection<some::vector<
                                some::scalar, some::indice_value<4>>>>,
      attribute<"cp_info", some::unbounded_collection<some::vector<
                               some::scalar, some::indice_value<25>>>>>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) positions(auto step)
    {
      auto& data = self()->data();
      using env_t = decltype(self()->env());
      //      using indice = typename env_t::indice;
      using scalar = typename env_t::scalar;

      auto& ids = storage::prop_values<system, "id">(data, step);
      auto& qs = storage::attr_values<system, "q">(data, step);

      attr<"pos_info">(*self()).clear();

      for (auto [id, q] : view::zip(ids, qs)) {
        attr<"pos_info">(*self()).push_back({(scalar)id, q[0], q[1], q[2]});
      }

      return algebra::matrix_view<algebra::unbounded_col_matrix<scalar, 4>>(
          attr<"pos_info">(*self()).data()->data(),
          attr<"pos_info">(*self()).size(),
          attr<"pos_info">(*self()).data()->size());
    }

    decltype(auto) velocities(auto step)
    {
      auto& data = self()->data();
      using env_t = decltype(self()->env());
      //      using indice = typename env_t::indice;
      using scalar = typename env_t::scalar;

      auto& ids = storage::prop_values<system, "id">(data, step);
      auto& velos = storage::attr_values<system, "velocity">(data, step);

      attr<"vel_info">(*self()).clear();

      for (auto [id, velo] : view::zip(ids, velos)) {
        attr<"vel_info">(*self()).push_back(
            {(scalar)id, velo[0], velo[1], velo[2]});
      }

      return algebra::matrix_view<algebra::unbounded_col_matrix<scalar, 4>>(
          attr<"vel_info">(*self()).data()->data(),
          attr<"vel_info">(*self()).size(),
          attr<"vel_info">(*self()).data()->size());
    }

    decltype(auto) contact_points(auto step)
    {
      auto& data = self()->data();
      using env_t = decltype(self()->env());
      using indice = typename env_t::indice;
      using scalar = typename env_t::scalar;

      /* /!\ first osi */
      auto& p0_v =
          storage::attr_values<osi, "p0_vector_assembled">(data, step)[0];

      auto& ys = storage::attr_values<interaction, "y">(data, step);
      auto& ydots = storage::attr_values<interaction, "ydot">(data, step);
      auto& lambdas = storage::attr_values<interaction, "lambda">(data, step);
      auto& nslaws = storage::attr_values<interaction, "nslaw">(data, step);
      auto& relations =
          storage::attr_values<interaction, "relation">(data, step);

      auto& ds1s = storage::prop_values<interaction, "ds1">(data, step);
      auto& ds2s = storage::prop_values<interaction, "ds2">(data, step);
      auto& activations =
          storage::prop_values<interaction, "activation">(data, step);

      attr<"cp_info">(*self()).clear();

      indice k = 0;
      for (auto [relation, nslaw, y, ydot, lambda, ds1, ds2, activation] :
           view::zip(relations, nslaws, ys, ydots, lambdas, ds1s, ds2s,
                     activations)) {
        auto hds1 = storage::handle(data, ds1);
        auto hds2 = storage::handle(data, ds2);

        using vect = std::decay_t<decltype(hds1.q())>; /* in 2D, 3 components:
                                                          translation 2 +
                                                          orientation 1 */

        if (activation) {
          auto index_ds1 =
              prop<"index">(hds1); /* cf one_step_intergator.hpp,
                                    * assemble_h_matrix_for_involved_ds => row
                                    * of p0_vector_assembled */
          auto p0 =
              algebra::get_vector(p0_v, index_ds1); /* in 2D, 2 components */

          vect c1 = {hds1.q()[0], hds1.q()[1], 0.};
          vect c2 = {hds2.q()[0], hds2.q()[1], 0.};

          vect cn;

          vect ca;

          vect cb;

          if (ds1 != ds2) {
            /* 2 disks */
            scalar dc2c1 = collision::distance(c2, c1);

            cn = (c2 - c1) / dc2c1;
            ca =
                c1 +
                storage::handle(data, storage::prop<"shape">(hds1)).radius() *
                    cn;
            cb =
                c2 -
                storage::handle(data, storage::prop<"shape">(hds2)).radius() *
                    cn;
          }
          else {
            variant::visit(
                data, relation,
                ground::overload(
                    /* disk / segment */
                    [&](storage::index<collision::disksegment_r, indice>
                            rrel) {
                      auto hrel = storage::handle(data, rrel);
                      /* cb is the proj point on the segment, computed a
                       * second time!
                       */
                      auto segment = storage::handle(data, hrel.segment());
                      const scalar t =
                          fmax(0, fmin(1, algebra::dot(c1 - segment.p1(),
                                                       segment.dp2p1()) /
                                              segment.length_sq()));
                      cb = segment.p1() + t * segment.dp2p1();
                      scalar dcbc1 = collision::distance(cb, c1);
                      cn = (cb - c1) / dcbc1;
                      ca = c1 + cn * storage::handle(
                                         data, storage::prop<"shape">(hds1))
                                         .radius();
                    },
                    /* disk / fixed disk */
                    [&](storage::index<collision::diskfdisk_r, indice> rrel) {
                      auto hrel = storage::handle(data, rrel);
                      auto tds = hrel.translated_disk_shape();

                      c2 = tds.translation();
                      scalar radius2 = tds.item().radius();

                      scalar dc2c1 = collision::distance(c2, c1);

                      cn = (c2 - c1) / dc2c1;
                      ca = c1 +
                           storage::handle(data, storage::prop<"shape">(hds1))
                                   .radius() *
                               cn;
                      cb = c2 - radius2 * cn;
                    },
                    [&](auto) {
                      throw(std::runtime_error(
                          "Error: not a disk-segment or disk-disk relation"));
                    }));
          };

          /* disk / segment */
          attr<"cp_info">(*self()).push_back(
              {storage::handle(data, nslaw).mu(),
               ca[0],
               ca[1],
               0., /* 2D */
               cb[0],
               cb[1],
               0., /* 2D */
               cn[0],
               cn[1],
               0, /* 2D */
               p0[0],
               p0[1],
               0, /* 2D */
               y[0],
               y[1],
               0, /* 2D */
               ydot[0],
               ydot[1],
               0, /* 2D */
               lambda[0],
               lambda[1],
               0, /* 2D */
               (scalar)k,
               (scalar)hds1.get(),
               (scalar)hds2.get()});
          k++;
        }
      }
      return algebra::matrix_view<algebra::unbounded_col_matrix<scalar, 25>>(
          attr<"cp_info">(*self()).data()->data(),
          attr<"cp_info">(*self()).size(),
          attr<"cp_info">(*self()).data()->size());
    }

    auto methods()
    {
      using env_t = decltype(self()->env());
      using indice = typename env_t::indice;

      return collect(
          method("positions", &interface<Handle>::positions<indice>),
          method("velocities", &interface<Handle>::velocities<indice>),
          method("contact_points",
                 &interface<Handle>::contact_points<indice>));
    }
  };
};
}  // namespace siconos::io
