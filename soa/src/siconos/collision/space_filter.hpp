#pragma once

#include <concepts>

#include "siconos/collision/collision_head.hpp"
#include "siconos/collision/diskdisk_r.hpp"
#include "siconos/collision/diskline_r.hpp"
#include "siconos/collision/point.hpp"
#include "siconos/collision/shape/disk.hpp"
#include "siconos/collision/shape/line.hpp"
#include "siconos/storage/pattern/base.hpp"
#include "siconos/storage/storage.hpp"
#include "siconos/utils/variant.hpp"

#define USE_DOUBLE
#include "CompactNSearch/CompactNSearch.h"

namespace siconos::collision {

/* the neighborhood engine */
template <typename... Points>
struct neighborhood
    : storage::data_holder<CompactNSearch::NeighborhoodSearch> {
  using points_t = gather<Points...>;

  using attributes = storage::pattern::cons_x<
      attribute<"point_set_id",
                some::array<some::indice,
                            some::indice_value<ground::size(points_t{})>>>,
      typename storage::data_holder<
          CompactNSearch::NeighborhoodSearch>::attributes>;

  template <typename Handle>
  struct interface : storage::data_holder<
                         CompactNSearch::NeighborhoodSearch>::
                         template interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) point_set_id()
    {
      return storage::attr<"point_set_id">(*self());
    };

    void create(auto radius)
    {
      self()->instance().reset(
          new CompactNSearch::NeighborhoodSearch(radius));
    }

    void add_point_sets(auto step)
    {
      auto& data = self()->data();

      using indice = typename decltype(self()->env())::indice;

      indice i = 0;

      auto& psid = storage::attr<"point_set_id">(*self());
      auto& instance = self()->instance();
      ground::for_each(points_t{}, [&data, &step, &i, &psid,
                                    &instance]<typename Point>(Point) {
        auto coords = storage::attr_values<Point, "coord">(data, step);
        psid[i++] =
            instance->add_point_set(coords.front().data(), coords.size());
      });
    }

    void update(auto step)
    {
      auto& data = self()->data();
      ground::for_each(points_t{}, [&data, &step]<typename Point>(Point) {
        for (auto point : storage::handles<Point>(data, step)) {
          point.update();
        }
      });
    }

    void search() { self()->instance()->find_neighbors(); };

    auto methods()
    {
      using env_t = decltype(self()->env());
      using indice = typename env_t::indice;
      using scalar = typename env_t::scalar;

      return collect(method("point_set_id", &interface<Handle>::point_set_id),
                     method("create", &interface<Handle>::create<scalar>),
                     method("add_point_sets",
                            &interface<Handle>::add_point_sets<indice>),
                     method("update", &interface<Handle>::update<indice>),
                     method("search", &interface<Handle>::search));
    }
  };
};

template <typename Topology, typename Neighborhood>
struct space_filter : item<> {
  using topology = Topology;
  using interaction = typename topology::interaction;
  using nslaw = typename interaction::nslaw;

  using neighborhood = Neighborhood;

  using attributes = gather<
      attribute<"topology", some::item_ref<topology>>,
      attribute<"neighborhood", some::item_ref<neighborhood>>,
      attribute<"nslaw", some::item_ref<nslaw>>,
      attribute<"diskdisk_r", some::item_ref<collision::diskdisk_r>>,
      attribute<"disklines",
                some::map<some::array<some::scalar, some::indice_value<3>>,
                          some::item_ref<diskline_r>>>>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) nslaw() { return storage::attr<"nslaw">(*self()); }
    decltype(auto) topology()
    {
      return storage::handle(self()->data(),
                             storage::attr<"topology">(*self()));
    }

    decltype(auto) neighborhood()
    {
      return storage::attr<"neighborhood">(*self());
    }

    decltype(auto) diskdisk_r()
    {
      return storage::attr<"diskdisk_r">(*self());
    }

    decltype(auto) disklines() { return storage::attr<"disklines">(*self()); }

    void make_points()
    {
      auto& data = self()->data();
      auto topo = self()->topology();

      using points_t = typename neighborhood::points_t;

      ground::for_each(points_t{}, [&]<typename Point>(Point) {
        using item_t = typename Point::item_t;
        if constexpr (std::derived_from<item_t, model::lagrangian_ds>) {
          auto all_ds = storage::handles<item_t>(data);
          for (auto ds : all_ds) {
            auto new_point = storage::add<Point>(data);
            new_point.item() = ds;
            new_point.update();

            auto dsgv = topo.dynamical_system_graphs()[0].add_vertex(ds);
            storage::prop<"vd">(ds) = dsgv;
          }
        }
        else if constexpr (std::derived_from<item_t,
                                             collision::shape::line>) {
          auto all_lines = storage::handles<item_t>(data);
          for (auto line : all_lines) {
            for (auto point_coord : line.points_coords()) {
              auto new_point = storage::add<Point>(data);
              new_point.item() = line;
              new_point.coord() = point_coord;
            }
          }
        }
      });
    }
    void update_index_set0()
    {
      using env = decltype(self()->env());

      auto& data = self()->data();
      auto topo = storage::handle(data, self()->topology());
      auto nslaw = self()->nslaw();
      auto diskdisk_r = self()->diskdisk_r();
      auto disklines = self()->disklines();

      auto ngbh = storage::handle(data, self()->neighborhood());
      using ngbh_t = typename std::decay_t<decltype(ngbh)>::type;
      using points_t = typename ngbh_t::points_t;

      auto& dsg0 = topo.dynamical_system_graphs()[0];
      auto& index_set0 = topo.interaction_graphs()[0];

      constexpr auto npointsets = ground::size(points_t{});
      ground::for_each(ground::range<npointsets>, [&](auto ip1) {
        auto p1 = points_t{}[ground::size_c<ip1>];
        using p1_t = decltype(p1);
        auto psid1 = ngbh.point_set_id()[ip1];

        ground::for_each(
            ground::range_c<std::size_t, ip1, npointsets>, [&](auto ip2) {
              auto p2 = points_t{}[ground::size_c<ip2>];
              using p2_t = decltype(p2);
              auto psid2 = ngbh.point_set_id()[ip2];

              auto& ps1 = ngbh.instance()->point_set(psid1);
              auto& ps2 = ngbh.instance()->point_set(psid2);

              for (size_t i = 0; i < ps1.n_points(); ++i) {
                auto pid1 = i;
                auto index_point1 = storage::index<p1_t, size_t>(pid1);
                auto handle_point1 = storage::handle(data, index_point1);
                auto body1 = storage::handle(data, handle_point1.item());

                for (size_t j = 0; j < ps2.n_neighbors(psid2, i); ++j) {
                  const unsigned int pid2 = ps2.neighbor(psid2, i, j);
                  auto index_point2 = storage::index<p2_t, size_t>(pid2);
                  auto handle_point2 = storage::handle(data, index_point2);
                  auto body2 = storage::handle(data, handle_point2.item());

                  using system1_t = typename p1_t::item_t;
                  using system2_t = typename p2_t::item_t;

                  if constexpr (std::derived_from<system1_t,
                                                  model::lagrangian_ds>) {
                    // proximity with another disk, only disks are dynamics
                    // check if interaction already exists

                    if constexpr (std::derived_from<system2_t,
                                                    model::lagrangian_ds>) {
                      auto& ds1 = body1;
                      auto& ds2 = body2;

                      auto& ds1d = storage::prop<"vd">(ds1);
                      // or dsg0.descriptor(ds1) with std::map
                      auto& ds2d =
                          storage::prop<"vd">(ds2);  // dsg0.descriptor(ds2);

                      // at most one edge between 2 ds !!
                      auto [edged, edge_exists] = dsg0.edge(ds1d, ds2d);

                      if (edge_exists) {
                        // keep this edge
                        auto inter =
                            storage::handle(data, dsg0.bundle(edged));
                        index_set0.color(storage::prop<"vd">(inter)) =
                            env::white_color;
                      }
                      else {
                        // create the edge
                        auto inter = topo.link(body1, body2);
                        inter.nslaw() = nslaw;  // one nslaw for the moment

                        index_set0.color(storage::prop<"vd">(inter)) =
                            env::white_color;
                        inter.relation() =
                            diskdisk_r;  // the diskdisk_r, need
                                         // only one relation!

                        // then flag edge (color ?) + remove all edges without
                        // the flag
                      }
                    }
                    else {
                      if constexpr (std::derived_from<system1_t,
                                                      model::lagrangian_ds>) {
                        // body2 is a static line
                        // for all self edges find the one with the
                        // corresponding line

                        auto line = storage::handle(data, body2);
                        auto a = line.a();
                        auto b = line.b();
                        auto c = line.c();

                        auto& ds1 = body1;
                        auto& ds1d = storage::prop<"vd">(ds1);

                        bool found = false;
                        for (auto [oei, oeiend] = dsg0.out_edges(ds1d);
                             (oei != oeiend); ++oei) {
                          if (dsg0.target(*oei) == ds1d) {
                            // self edge
                            // is it this line ?
                            auto inter =
                                storage::handle(data, dsg0.bundle(*oei));

                            auto vds =
                                storage::prop_values<interaction, "vd">(data,
                                                                        0);

                            if (siconos::variant::visit(
                                    data, inter.relation(),
                                    ground::overload(
                                        [&]<match::handle<diskline_r>
                                                DiskLineR>(DiskLineR rel) {
                                          auto line = storage::handle(
                                              data, rel.line());
                                          return line.a() == a &&
                                                 line.b() == b &&
                                                 line.c() == c;
                                        },
                                        []<bool flag = false>(auto) {
                                          assert(flag);
                                          return false;
                                          // static_assert(flag,
                                          //               "should not
                                          //               happen");
                                        }))) {
                              // keep this interaction
                              index_set0.color(storage::prop<"vd">(inter)) =
                                  env::white_color;
                              found = true;  // found the edge
                            }
                            else {  // found an edge corresponding to another
                                    // line, nothing to do
                            }
                          }
                          else {  // not a self edge, nothing to do
                          }
                        }

                        auto vds =
                            storage::prop_values<interaction, "vd">(data, 0);

                        if (!found) {
                          // create the edge
                          auto inter = topo.link(body1);
                          inter.nslaw() = nslaw;  // one nslaw for the moment

                          if (auto search = disklines.find({a, b, c});
                              search != disklines.end()) {
                            inter.relation() = search->second;
                          }
                          else {
                            auto dl = storage::add<diskline_r>(data);
                            auto line = storage::handle(data, dl.line());
                            line.a() = a;
                            line.b() = b;
                            line.c() = c;

                            inter.relation() = dl;
                          }
                          index_set0.color(storage::prop<"vd">(inter)) =
                              env::white_color;
                        }
                        else {
                          break;
                        }
                      }
                    }
                  }
                }
              }
            });
      });
      auto hh = storage::handles<interaction>(data);
      for (auto inter : hh) {
        auto& inter_vd = storage::prop<"vd">(inter);
        if (!inter_vd) {
          // fix: this case should not exist
          storage::remove(data, inter);
        }
        else if (index_set0.color(inter_vd) != env::white_color) {
          // remove interaction
          auto& ds1d = storage::prop<"ds1d">(inter);
          //          auto& ds2d = storage::prop<"ds2d">(inter);

          // interaction graph
          index_set0.remove_vertex(inter_vd);

          // dynamical system graph
          bool done = false;
          for (auto [oei, oeiend] = dsg0.out_edges(ds1d);
               !done && (oei != oeiend); ++oei) {
            if (dsg0.bundle(*oei).get() == inter.get()) {
              // only one edge for an interaction
              dsg0.remove_edge(*oei);
              // iterator invalidated, we must escape
              done = true;
            }
          }

          // with move_back
          storage::remove(data, inter);
        }
      }
    }

    auto methods()
    {
      //      using env_t = decltype(self()->env());
      //      using indice = typename env_t::indice;
      //      using scalar = typename env_t::scalar;

      return collect(
          method("make_points", &interface<Handle>::make_points),
          method("update_index_set0", &interface<Handle>::update_index_set0));
    }
  };
};
}  // namespace siconos::collision
