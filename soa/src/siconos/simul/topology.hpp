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

  using attributes =
      gather<attribute<"dynamical_system_graphs",
                       some::bounded_collection<
                           some::graph<some::item_ref<dynamical_system>,
                                       some::item_ref<interaction>>,
                           std::integral_constant<int, 1>>>,
             attribute<"interaction_graphs",
                       some::bounded_collection<
                           some::graph<some::item_ref<interaction>,
                                       some::item_ref<dynamical_system>>,
                           std::integral_constant<int, 2>>>>;

  using properties = gather<
      storage::attached<dynamical_system, symbol<"involved">, some::boolean>,
      storage::attached<dynamical_system, symbol<"index">, some::indice>,
      storage::attached<dynamical_system, symbol<"id">, some::indice>,
      storage::attached<dynamical_system, symbol<"p0">,
                        some::array<some::vector<some::scalar, dof>,
                                    std::integral_constant<int, 2>>>,
      storage::attached<interaction, symbol<"ds1d">,
                        some::vdescriptor<typename get_attr_t<
                            attributes, "dynamical_system_graphs">::type>>,
      storage::attached<interaction, symbol<"ds2d">,
                        some::vdescriptor<typename get_attr_t<
                            attributes, "dynamical_system_graphs">::type>>,
      storage::attached<interaction, symbol<"nds">, some::indice>,
      storage::attached<interaction, symbol<"ds1">,
                        some::item_ref<dynamical_system>>,
      storage::attached<interaction, symbol<"ds2">,
                        some::item_ref<dynamical_system>>,
      storage::attached<interaction, symbol<"activation">, some::boolean>,
      storage::attached<dynamical_system, symbol<"vd">,
                        some::vdescriptor<typename get_attr_t<
                            attributes, "dynamical_system_graphs">::type>>,
      storage::attached<interaction, symbol<"vd">,
                        some::vdescriptor<typename get_attr_t<
                            attributes, "interaction_graphs">::type>>>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    decltype(auto) dynamical_system_graphs()
    {
      return attr<"dynamical_system_graphs">(*self());
    };
    decltype(auto) interaction_graphs()
    {
      return attr<"interaction_graphs">(*self());
    };

    using default_interface<Handle>::self;

    template <match::handle<dynamical_system> Hds>
    decltype(auto) link(Hds ds)
    {
      auto &data = self()->data();

      auto &dsg0 = self()->dynamical_system_graphs()[0];
      auto &ig0 = self()->interaction_graphs()[0];
      auto inter = storage::add<interaction>(data);
      auto dsgv = dsg0.add_vertex(ds);
      auto [dsg0_ed, ig0_vd] = dsg0.add_edge(dsgv, dsgv, inter, ig0);

      attr<"y">(inter).setZero();
      attr<"ydot">(inter).setZero();

      prop<"vd">(ds) = dsgv;
      prop<"vd">(inter) = ig0_vd;
      prop<"nds">(inter) = 1;
      prop<"ds1">(inter) = ds;
      prop<"ds2">(inter) = ds;

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

      auto &dsg0 = self()->dynamical_system_graphs()[0];
      auto &ig0 = self()->interaction_graphs()[0];

      auto inter = storage::add<interaction>(data);
      auto dsgv1 = dsg0.add_vertex(ds1);
      auto dsgv2 = dsg0.add_vertex(ds2);
      auto [dsg0_ed, ig0_vd] = dsg0.add_edge(dsgv1, dsgv2, inter, ig0);
      // dsg0_ed may be invalidated on edge removal so it is useless to store
      // it

      attr<"y">(inter).setZero();
      attr<"ydot">(inter).setZero();

      ds1.property(symbol<"vd">{}) = dsgv1;
      ds2.property(symbol<"vd">{}) = dsgv2;

      inter.property(symbol<"nds">{}) = 2;
      inter.property(symbol<"vd">{}) = ig0_vd;
      inter.property(symbol<"ds1">{}) = ds1;
      inter.property(symbol<"ds2">{}) = ds2;

      inter.property(symbol<"ds1d">{}) = dsgv1;
      inter.property(symbol<"ds2d">{}) = dsgv2;

      return inter;
    };
  };
};
}  // namespace siconos::simul
