#pragma once

#include "siconos/simul/simul_head.hpp"
#include "siconos/utils/print.hpp"
#include "siconos/utils/range.hpp"

namespace siconos::simul {

template <typename DynamicalSystem, typename Interaction>
struct topology : item<> {
  using dof = some::indice_parameter<"dof">;

  using dynamical_system = DynamicalSystem;
  using interaction = Interaction;

  using nslaw = typename interaction::nslaw;
  using nslaw_size = some::indice_value<nslaw::size>;

  using attributes = gather<
      attribute<"system_id",
                some::map<some::indice, some::item_ref<dynamical_system>>>>;

  using properties = gather<
      storage::attached<dynamical_system, symbol<"involved">, some::boolean>,
      storage::attached<dynamical_system, symbol<"index">, some::indice>,
      storage::attached<dynamical_system, symbol<"id">, some::indice>,
      storage::attached<dynamical_system, symbol<"p0">,
                        some::array<some::vector<some::scalar, dof>,
                                    std::integral_constant<int, 2>>>,

      storage::attached<interaction, symbol<"nds">, some::indice>,
      storage::attached<interaction, symbol<"ds1">,
                        some::item_ref<dynamical_system>>,
      storage::attached<interaction, symbol<"ds2">,
                        some::item_ref<dynamical_system>>,
      storage::attached<interaction, symbol<"ydot_backup">,
                        some::vector<some::scalar, nslaw_size>>,
      storage::attached<interaction, symbol<"activation">, some::boolean>>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    template <match::handle<dynamical_system> Hds>
    decltype(auto) link(Hds ds)
    {
      auto &data = self()->data();

      auto inter = storage::add<interaction>(data);
      //      auto dsgv = dsg0.add_vertex(ds);
      //      auto [dsg0_ed, ig0_vd] = dsg0.add_edge(dsgv, dsgv, inter, ig0);

      attr<"y">(inter).setZero();
      attr<"ydot">(inter).setZero();
      attr<"lambda">(inter).setZero();

      //      prop<"vd">(ds) = dsgv;
      //      prop<"vd">(inter) = ig0_vd;
      prop<"nds">(inter) = 1;
      prop<"ds1">(inter) = ds;
      prop<"ds2">(inter) = ds;

      //      prop<"ds1d">(inter) = dsgv;
      //      prop<"ds2d">(inter) = dsgv;
      //      dsg0.update_vertices_indices();
      //      dsg0.update_edges_indices();

      //      ig0.update_vertices_indices();
      //      ig0.update_edges_indices();

      /*      print("dsg0:\n:");
            print("=========\n");
            dsg0.display();
            print("=========\n");

            print("ig0:\n");
            print("=========\n");
            ig0.display();
            print("=========\n");*/

      return inter;
    };

    template <match::handle<dynamical_system> Hds>
    decltype(auto) link(Hds ds1, Hds ds2)
    {
      auto &data = self()->data();

      auto inter = storage::add<interaction>(data);
      //      auto dsgv1 = dsg0.add_vertex(ds1);
      //      auto dsgv2 = dsg0.add_vertex(ds2);
      //      auto [dsg0_ed, ig0_vd] = dsg0.add_edge(dsgv1, dsgv2, inter,
      //      ig0);
      // dsg0_ed may be invalidated on edge removal so it is useless to store
      // it

      attr<"y">(inter).setZero();
      attr<"ydot">(inter).setZero();
      attr<"lambda">(inter).setZero();

      //      ds1.property(symbol<"vd">{}) = dsgv1;
      //      ds2.property(symbol<"vd">{}) = dsgv2;

      inter.property(symbol<"nds">{}) = 2;
      //      inter.property(symbol<"vd">{}) = ig0_vd;
      inter.property(symbol<"ds1">{}) = ds1;
      inter.property(symbol<"ds2">{}) = ds2;

      //      inter.property(symbol<"ds1d">{}) = dsgv1;
      //      inter.property(symbol<"ds2d">{}) = dsgv2;

      return inter;
    };

    void set_dynamical_system_id(auto hds, auto id)
    {
      attr<"system_id">(*self())[id] = hds.index_cast();
    }

    decltype(auto) dynamical_system(auto id)
    {
      return storage::handle(self()->data(), attr<"system_id">(*self())[id]);
    }

    auto methods()
    {
      using data_t = std::decay_t<decltype(self()->data())>;
      using env_t = decltype(self()->env());
      using indice = typename env_t::indice;
      using ds_t = DynamicalSystem;
      using hds_t = storage::handle<ds_t, indice, data_t>;

      return collect(
          method("set_dynamical_system_id",
                 &interface<Handle>::set_dynamical_system_id<hds_t, indice>),
          method("dynamical_system",
                 &interface<Handle>::dynamical_system<indice>));
    }
  };
};
}  // namespace siconos::simul
