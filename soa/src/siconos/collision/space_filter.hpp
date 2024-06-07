#pragma once

#include <concepts>

#include "siconos/collision/collision_head.hpp"
#include "siconos/collision/diskdisk_r.hpp"
#include "siconos/collision/diskfdisk_r.hpp"
#include "siconos/collision/disksegment_r.hpp"
#include "siconos/collision/point.hpp"
#include "siconos/collision/shape/disk.hpp"
#include "siconos/collision/shape/segment.hpp"
#include "siconos/storage/pattern/base.hpp"
#include "siconos/storage/storage.hpp"
#include "siconos/utils/print.hpp"
#include "siconos/utils/variant.hpp"

#define USE_DOUBLE

#ifdef WITH_GPU
#include "cuNSearch/cuNSearch.h"
namespace CompactNSearch = cuNSearch;
#else
#include "CompactNSearch/CompactNSearch.h"
#endif

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
      using env_t = decltype(self()->env());
      using indice = typename env_t::indice;

      auto& data = self()->data();

      indice i = 0;

      auto& psid = storage::attr<"point_set_id">(*self());
      auto& instance = self()->instance();
      ground::for_each(points_t{}, [&data, &step, &i, &psid,
                                    &instance]<typename Point>(Point) {
        auto& coords = storage::attr_values<Point, "coord">(data, step);

        psid[i++] =
            instance->add_point_set(coords.front().data(), coords.size());
      });
    }

    void set_active(auto ps1_id, auto ps2_id, auto value)
    {
      self()->instance()->set_active((unsigned int)ps1_id,
                                     (unsigned int)ps2_id, value);
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

    void sort()
    {
      using env_t = decltype(self()->env());
      using indice = typename env_t::indice;

      auto& data = self()->data();
      auto& instance = self()->instance();

      instance->z_sort();

      indice i = 0;
      ground::for_each(
          points_t{}, [&instance, &data, &i]<typename Point>(Point p) {
            auto ps = instance->point_set(i++);
            // apply function only if some points exist
            if (ps.n_points() > 0) {
              storage::apply_fun(data, p, [&ps]<typename Array>(Array& a) {
                  ps.sort_field(a.data());
                });
            }
          });
    }
    auto methods()
    {
      using env_t = decltype(self()->env());

      using indice = typename env_t::indice;
      using scalar = typename env_t::scalar;

      return collect(
          method("point_set_id", &interface<Handle>::point_set_id),
          method("create", &interface<Handle>::create<scalar>),
          method("add_point_sets",
                 &interface<Handle>::add_point_sets<indice>),
          method("update", &interface<Handle>::update<indice>),
          method("set_active",
                 &interface<Handle>::set_active<indice, indice, bool>),
          method("search", &interface<Handle>::search),
          method("sort", &interface<Handle>::sort));
    }
  };
};

template <typename Topology, typename Neighborhood>
struct space_filter : item<> {
  using topology = Topology;
  using dynamical_system = typename topology::dynamical_system;
  using interaction = typename topology::interaction;
  using nslaw = typename interaction::nslaw;

  using neighborhood = Neighborhood;

  using attributes = gather<
      attribute<"topology", some::item_ref<topology>>,
      attribute<"neighborhood", some::item_ref<neighborhood>>,
      attribute<"nslaw", some::item_ref<nslaw>>,
      attribute<"diskdisk_r", some::item_ref<diskdisk_r>>,
      attribute<"disksegments",
                some::map<some::array<some::scalar, some::indice_value<4>>,
                          some::item_ref<disksegment_r>>>,
      attribute<"diskfdisks",
                some::map<some::array<some::scalar, some::indice_value<2>>,
                          some::item_ref<diskfdisk_r>>>>;

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
      return storage::handle(self()->data(),
                             storage::attr<"neighborhood">(*self()));
    }

    decltype(auto) diskdisk_r()
    {
      return storage::attr<"diskdisk_r">(*self());
    }

    decltype(auto) disksegments()
    {
      return storage::attr<"disksegments">(*self());
    }

    decltype(auto) diskfdisks()
    {
      return storage::attr<"diskfdisks">(*self());
    }

    void remove_static_item(auto step, auto item_handle)
    {
      using env_t = decltype(self()->env());
      using indice = typename env_t::indice;
      using item_t = typename std::decay_t<decltype(item_handle)>::type;
      using point_t = collision::point<item_t>;
      using points_t = typename neighborhood::points_t;

      auto& data = self()->data();

      auto& points_flags = storage::attr_values<point_t, "flag">(data, step);
      auto& points_items = storage::attr_values<point_t, "item">(data, step);
      auto& points_coords =
          storage::attr_values<point_t, "coord">(data, step);

      auto ps_indx = ground::index_of<point_t>(ground::std_tuple(points_t{}));

      // ground::index_if(
      //   points_t{},
      //   ground::equal.to(collision::point<item_t>));

      // first remove item
      storage::remove(data, item_handle);

      // find all associated points
      for (auto [flag] : view::zip(points_flags)) {
        flag = false;
      }

      for (auto [flag, pitem] : view::zip(points_flags, points_items)) {
        if (item_handle.get() == pitem.get()) {
          flag = true;
        };
      }

      // remove points
      auto ff = std::ranges::find(points_flags, true);

      while (ff != points_flags.end()) {
        auto ff_index = ff - points_flags.begin();
        auto point =
            storage::handle(data, storage::index<point_t, indice>(ff_index));
        storage::remove(data, point);

        auto remaining_points_flags =
            std::ranges::subrange(ff, points_flags.end());
        auto rff = std::ranges::find(remaining_points_flags, true);

        ff = ff + (rff - remaining_points_flags.begin());
      };

      neighborhood().instance()->resize_point_set(
          ps_indx, points_coords.front().data(), points_coords.size());
    }

    void insert_disksegment_r(auto dl)
    {
      auto hdl = storage::handle(self()->data(), dl);
      storage::attr<"disksegments">(
          *self())[{hdl.segment().x1(), hdl.segment().x2(),
                    hdl.segment().y1(), hdl.segment().y2()}] = dl;

      // for (auto hds : storage::handles<dynamical_system>(self()->data())) {
      //   auto inter = self()->topology().link(hds);
      //   inter.relation() = dl;
      //   inter.nslaw() = nslaw();
      // }
    }

    void insert_diskfdisk_r(auto tds)
    {
      auto htds =
          storage::handle(self()->data(), tds).translated_disk_shape();
      storage::attr<"diskfdisks">(
          *self())[{htds.translation()[0], htds.translation()[1]}] = tds;

      // for (auto hds : storage::handles<dynamical_system>(self()->data())) {
      //   auto inter = self()->topology().link(hds);
      //   inter.relation() = dl;
      //   inter.nslaw() = nslaw();
      // }
    }

    void make_points()
    {
      auto& data = self()->data();

      using points_t = typename neighborhood::points_t;

      ground::for_each(points_t{}, [&]<typename Point>(Point) {
        using item_t = typename Point::item_t;
        if constexpr (std::derived_from<item_t, model::lagrangian_ds>) {
          auto all_ds = storage::handles<item_t>(data);
          for (auto ds : all_ds) {
            //            print("add disk point : {}\n", ds.get());
            auto new_point = storage::add<Point>(data);
            new_point.item() = ds;
            new_point.update();
          }
        }
        else if constexpr (std::derived_from<item_t,
                                             collision::shape::segment>) {
          auto all_segments = storage::handles<item_t>(data);
          for (auto segment : all_segments) {
            for (auto point_coord : segment.points_coords()) {
              //              print("segment point for {},{},{},{} : {},{}\n",
              //              segment.x1(),
              //                    segment.x2(), segment.y1(), segment.y2(),
              //                    point_coord(0), point_coord(1));
              auto new_point = storage::add<Point>(data);
              new_point.item() = segment;
              new_point.coord() = point_coord;
            }
          }
        }
        else if constexpr (std::derived_from<item_t,
                                             collision::translated<
                                                 collision::shape::disk>>) {
          auto all_fdisks = storage::handles<item_t>(data);
          for (auto fdisk : all_fdisks) {
            auto new_point = storage::add<Point>(data);
            new_point.item() = fdisk;
            new_point.coord() = fdisk.translation();
          }
        };
      });
    }

    decltype(auto) make_ipair(auto ds1, auto ds2)
    {
      auto i1 = ds1.get();
      auto i2 = ds2.get();
      if (i1 < i2) {
        return ground::make_pair(i1, i2);
      }
      else {
        return ground::make_pair(i2, i1);
      }
    }
    void update_index_set0(auto step)
    {
      using env = decltype(self()->env());
      using indice = typename env::indice;
      using scalar = typename env::scalar;

      auto& data = self()->data();
      auto topo = storage::handle(data, self()->topology());
      auto nslaw = self()->nslaw();
      auto diskdisk_r = self()->diskdisk_r();
      auto disksegments = self()->disksegments();
      auto diskfdisks = self()->diskfdisks();

      auto ngbh = storage::handle(data, self()->neighborhood());
      using ngbh_t = typename std::decay_t<decltype(ngbh)>::type;
      using points_t = typename ngbh_t::points_t;

      auto ds_ds_prox = std::map<ground::pair<indice, indice>,
                                 storage::index<interaction, indice>>();

      auto ds_segment_prox =
          std::map<ground::pair<indice, std::array<scalar, 4>>,
                   storage::index<interaction, indice>>();

      auto ds_fdisk_prox =
          std::map<ground::pair<indice, std::array<scalar, 2>>,
                   storage::index<interaction, indice>>();

      auto& ds1s = storage::prop_values<interaction, "ds1">(data, step);
      auto& ds2s = storage::prop_values<interaction, "ds2">(data, step);
      auto interactions = storage::handles<interaction>(data, step);

      // build (ds ds -> inter) & (ds segment -> inter) maps
      for (auto [ds1, ds2, inter] : view::zip(ds1s, ds2s, interactions)) {
        if (ds1 != ds2) {
          ds_ds_prox[make_ipair(ds1, ds2)] = inter;
        }
        else {
          siconos::variant::visit(
              data, inter.relation(),
              ground::overload(
                  // https://stackoverflow.com/questions/46114214/lambda-implicit-capture-fails-with-variable-declared-from-structured-binding
                  // capture by value ok with handles
                  [&, ds1 = ds1,
                   inter = inter]<match::handle<disksegment_r> DiskSegmentR>(
                      DiskSegmentR rel) {
                    auto segment = rel.segment();
                    auto coefs = std::array{segment.x1(), segment.x2(),
                                            segment.y1(), segment.y2()};
                    ds_segment_prox[ground::make_pair(ds1.get(), coefs)] =
                        inter;
                  },
                  [&, ds1 = ds1,
                   inter = inter]<match::handle<diskfdisk_r> DiskFdiskR>(
                      DiskFdiskR rel) {
                    auto fdisk = rel.translated_disk_shape();
                    auto& translat = fdisk.translation();
                    auto coefs = std::array{translat[0], translat[1]};
                    ds_fdisk_prox[ground::make_pair(ds1.get(), coefs)] =
                        inter;
                  },
                  []<bool flag = false>(auto) {
                    assert(flag);
                    // should send an exception here
                    // static_assert(flag,
                    //               "should not
                    //               happen");
                  }));
        }
      }

      auto& activations =
          storage::prop_values<interaction, "activation">(data, 0);

      int activations_size = activations.size();
      if (activations_size > 0) {
        for (auto [activation] : view::zip(activations)) {
          activation = false;
        }
      }

      constexpr auto npointsets = ground::size(points_t{});
      ground::for_each(
          ground::range<npointsets - ground::size_c<1_c>>, [&](auto ip1) {
            auto p1 = points_t{}[ground::size_c<ip1>];
            using p1_t = decltype(p1);
            auto psid1 = ngbh.point_set_id()[ip1];

            ground::for_each(
                ground::range_c<std::size_t, ip1, npointsets>, [&](auto ip2) {
                  auto p2 = points_t{}[ground::size_c<ip2>];
                  using p2_t = decltype(p2);
                  auto psid2 = ngbh.point_set_id()[ip2];

                  auto& ps1 = ngbh.instance()->point_set(psid1);
                  //              auto& ps2 =
                  //              ngbh.instance()->point_set(psid2);

                  for (size_t i = 0; i < ps1.n_points(); ++i) {
                    auto pid1 = i;
                    auto index_point1 = storage::index<p1_t, size_t>(pid1);
                    auto handle_point1 = storage::handle(data, index_point1);
                    auto body1 = storage::handle(data, handle_point1.item());

                    for (size_t j = 0; j < ps1.n_neighbors(psid2, i); ++j) {
                      const unsigned int pid2 = ps1.neighbor(psid2, i, j);

                      // print("pid2 : {}\n", pid2);
                      auto index_point2 = storage::index<p2_t, size_t>(pid2);
                      auto handle_point2 =
                          storage::handle(data, index_point2);
                      auto body2 =
                          storage::handle(data, handle_point2.item());
                      using system1_t = typename p1_t::item_t;
                      using system2_t = typename p2_t::item_t;

                      // print("point1 {},{},{}\n",
                      //       ground::type_name<system1_t>().c_str(),
                      //       handle_point1.coord()(0),
                      //       handle_point1.coord()(1));
                      // print("point2 {},{},{}\n",
                      //       ground::type_name<system2_t>().c_str(),
                      //       handle_point2.coord()(0),
                      //       handle_point2.coord()(1));

                      if constexpr (std::derived_from<system1_t,
                                                      model::lagrangian_ds>) {
                        // proximity with another disk, only disks are
                        // dynamics check if interaction already exists

                        if constexpr (std::derived_from<
                                          system2_t, model::lagrangian_ds>) {
                          auto& ds1 = body1;
                          auto& ds2 = body2;

                          // at most one edge between 2 ds !!
                          auto find_inter =
                              ds_ds_prox.find(make_ipair(ds1, ds2));
                          if (find_inter != ds_ds_prox.end()) {
                            // keep this edge
                            auto inter = storage::handle(
                                data, std::get<1>(*find_inter));
                            storage::prop<"activation">(inter) = true;
                          }
                          else {
                            // create the edge
                            auto inter = topo.link(body1, body2);
                            inter.nslaw() =
                                nslaw;  // one nslaw for the moment

                            storage::prop<"activation">(inter) = true;

                            inter.relation() =
                                diskdisk_r;  // the diskdisk_r, need
                                             // only one relation!
                            ds_ds_prox[make_ipair(ds1, ds2)] = inter;
                          }
                        }
                        else {
                          if constexpr (std::derived_from<
                                            system1_t,
                                            model::lagrangian_ds>) {
                            if constexpr (std::derived_from<
                                              system2_t,
                                              collision::shape::segment>) {
                              // body2 is a static segment
                              // for all self edges find the one with the
                              // corresponding segment

                              auto segment = storage::handle(data, body2);
                              auto x1 = segment.x1();
                              auto x2 = segment.x2();
                              auto y1 = segment.y1();
                              auto y2 = segment.y2();

                              // print("PROXIMITY with {},{},{}\n", x1, x2,
                              // y1, y2);

                              auto& ds1 = body1;

                              auto find_inter =
                                  ds_segment_prox.find(ground::make_pair(
                                      ds1.get(), std::array{x1, x2, y1, y2}));

                              if (find_inter != ds_segment_prox.end()) {
                                auto inter = storage::handle(
                                    data, std::get<1>(*find_inter));
                                // keep this interaction
                                storage::prop<"activation">(inter) = true;

                                // print("interaction FOUND for {},{},{}\n",
                                // a, b,
                                //       c);
                              }
                              else {
                                // print("interaction NOT FOUND for
                                // {},{},{}\n", a,
                                //       b, c);
                                // create the edge
                                auto inter = topo.link(body1);
                                inter.nslaw() =
                                    nslaw;  // one nslaw for the moment

                                if (auto search =
                                        disksegments.find({x1, x2, y1, y2});
                                    search != disksegments.end()) {
                                  inter.relation() = search->second;
                                }
                                else {
                                  auto dl = storage::add<disksegment_r>(data);
                                  auto segment =
                                      storage::handle(data, dl.segment());
                                  segment.x1() = x1;
                                  segment.x2() = x2;
                                  segment.y1() = y1;
                                  segment.y2() = y2;

                                  inter.relation() = dl;
                                  disksegments[{x1, x2, y1, y2}] = dl;
                                }
                                storage::prop<"activation">(inter) = true;
                                ds_segment_prox[ground::make_pair(
                                    ds1.get(), std::array{x1, x2, y1, y2})] =
                                    inter;
                              }
                            }
                            else if constexpr (
                                std::derived_from<
                                    system2_t, collision::translated<
                                                   collision::shape::disk>>) {
                              auto fdisk = storage::handle(data, body2);
                              auto& translat = fdisk.translation();
                              auto coefs =
                                  std::array{translat[0], translat[1]};

                              auto& ds1 = body1;

                              auto find_inter = ds_fdisk_prox.find(
                                  ground::make_pair(ds1.get(), coefs));

                              if (find_inter != ds_fdisk_prox.end()) {
                                auto inter = storage::handle(
                                    data, std::get<1>(*find_inter));
                                storage::prop<"activation">(inter) = true;
                              }
                              else {
                                auto inter = topo.link(body1);
                                inter.nslaw() = nslaw;

                                if (auto search = diskfdisks.find(coefs);
                                    search != diskfdisks.end()) {
                                  inter.relation() = search->second;
                                }
                                else {
                                  auto dfd = storage::add<diskfdisk_r>(data);
                                  auto tds = storage::handle(
                                      data, dfd.translated_disk_shape());
                                  tds.translation() = translat;

                                  inter.relation() = dfd;
                                  diskfdisks[coefs] = dfd;
                                }
                                storage::prop<"activation">(inter) = true;
                                ds_fdisk_prox[ground::make_pair(
                                    ds1.get(), coefs)] = inter;
                              }
                            }
                          }
                        }
                      }
                    }
                  };
                });
          });

      //      print("BEFORE REMOVAL: size of indexset0: {}\n",
      //      activations.size());
      for (auto [activation, inter] : view::zip(activations, interactions)) {
        if (!activation) {
          //          print("START REMOVE interaction {}\n", inter.get());

          if (storage::prop<"ds1">(inter) != storage::prop<"ds2">(inter)) {
            auto finter = ds_ds_prox.find(make_ipair(
                storage::prop<"ds1">(inter), storage::prop<"ds2">(inter)));

            assert(inter.get() == std::get<1>(*finter).get());

            ds_ds_prox.erase(finter);
            //            print("  REMOVE ds ds interaction between {}
            //            {}\n",
            //                  storage::prop<"ds1">(inter).get(),
            //                  storage::prop<"ds2">(inter).get());
          }
          else {
            siconos::variant::visit(
                data, inter.relation(),
                ground::overload(
                    [&, inter =
                            inter]<match::handle<disksegment_r> DiskSegmentR>(
                        DiskSegmentR rel) {
                      auto segment = rel.segment();
                      auto coefs = std::array{segment.x1(), segment.x2(),
                                              segment.y1(), segment.y2()};
                      auto finter = ds_segment_prox.find(ground::make_pair(
                          storage::prop<"ds1">(inter).get(), coefs));
                      ds_segment_prox.erase(finter);
                    },
                    [&, inter = inter]<match::handle<diskfdisk_r> DiskFdiskR>(
                        DiskFdiskR rel) {
                      auto& translat =
                          rel.translated_disk_shape().translation();
                      auto coefs = std::array{translat[0], translat[1]};
                      auto finter = ds_fdisk_prox.find(ground::make_pair(
                          storage::prop<"ds1">(inter).get(), coefs));
                      ds_fdisk_prox.erase(finter);
                    },
                    []<bool flag = false>(auto) {
                      assert(flag);
                      // should send an exception here
                      // static_assert(flag,
                      //               "should not
                      //               happen");
                    }));
            //            print("  REMOVE ds segment interaction between {}
            //            {},{},{})\n",
            //                  storage::prop<"ds1">(inter).get(),
            //                  segmentcoefs[0], segmentcoefs[1],
            //                  segmentcoefs[2], segmentcoefs[3]);

            ;
          }
        }
      }

      auto fact = std::ranges::find(activations, false);

      //      print("  START REMOVE interactions\n");

      while (fact != activations.end()) {
        auto fact_index = fact - activations.begin();
        //        print("  activation of {} is false\n", fact_index);
        auto inter = storage::handle(
            data, storage::index<interaction, indice>(fact_index));

        // with move_back : order is modified

        //        print("  remove interaction {}\n", inter.get());
        storage::remove(data, inter);

        // XXX remove from ds_ds_prox or ds_segment_prox

        // activations has been modified, search first false element
        // starting at current position
        auto remaining_activations =
            std::ranges::subrange(fact, activations.end());

        auto ifact = std::ranges::find(remaining_activations, false);

        fact = fact + (ifact - remaining_activations.begin());
        // fact = std::ranges::find(activations, false);

        //        print("  find new false at : {}\n", fact -
        //        activations.begin());
      }

      //      print("AFTER REMOVAL: size of indexset0: {}\n",
      //      activations.size()); print("END of interactions removal\n");
      //      print("size of ds ds map: {}\n", ds_ds_prox.size());
      //      print("size of ds segment map: {}\n", ds_segment_prox.size());
    }

    auto methods()
    {
      auto& data = self()->data();
      using env_t = decltype(self()->env());
      using indice = typename env_t::indice;
      // using scalar = typename env_t::scalar;
      using disksegment_r_t =
          storage::index<collision::disksegment_r, indice>;
      using diskfdisk_r_t = storage::index<collision::diskfdisk_r, indice>;
      using segment_handle_t = std::decay_t<decltype(storage::handle(
          data, storage::index<collision::shape::segment, indice>{}))>;

      return collect(
          method("make_points", &interface<Handle>::make_points),
          method("update_index_set0",
                 &interface<Handle>::update_index_set0<indice>),
          method("insert_disksegment_r",
                 &interface<Handle>::insert_disksegment_r<disksegment_r_t>),
          method("insert_diskfdisk_r",
                 &interface<Handle>::insert_diskfdisk_r<diskfdisk_r_t>),
          method("remove_static_segment",
                 &interface<Handle>::template remove_static_item<
                     indice, segment_handle_t>));
    }
  };
};
}  // namespace siconos::collision
