#pragma once

#include "siconos/siconos.hpp"
#include "siconos/utils/pattern.hpp"

namespace siconos {

template <typename DynamicalSystem, typename Interaction>
struct topology : item<> {
  using dof = some::indice_parameter<"dof">;

  using dynamical_system = DynamicalSystem;
  using interaction = Interaction;

  struct dynamical_system_graphs
      : some::bounded_collection<some::graph<some::item_ref<dynamical_system>,
                                             some::item_ref<interaction> >,
                                 std::integral_constant<int, 1> >,
        access<dynamical_system_graphs> {};

  struct interaction_graphs
      : some::bounded_collection<
            some::graph<some::item_ref<interaction>,
                        some::item_ref<dynamical_system> >,
            std::integral_constant<int, 2> >,
        access<interaction_graphs> {};

  // number of involved ds
  struct nids : some::indice, access<nids> {};

  using attributes =
      gather<dynamical_system_graphs, interaction_graphs, nids>;

  using properties = gather<
      attached_storage<dynamical_system, symbol<"involved">, some::boolean>,
      attached_storage<dynamical_system, symbol<"index">, some::indice>,
      attached_storage<dynamical_system, symbol<"p0">,
                       some::vector<some::vector<some::scalar, dof>,
                                    std::integral_constant<int, 2> > >,
      attached_storage<dynamical_system, symbol<"velocity">,
                       some::vector<some::vector<some::scalar, dof>,
                                    std::integral_constant<int, 2> > >,
      attached_storage<dynamical_system, symbol<"q">,
                       some::vector<some::vector<some::scalar, dof>,
                                    std::integral_constant<int, 2> > >,
      attached_storage<interaction, symbol<"nds">, some::indice>,
      attached_storage<interaction, symbol<"ds1">,
                       some::item_ref<dynamical_system> >,
      attached_storage<interaction, symbol<"ds2">,
                       some::item_ref<dynamical_system> >,
      attached_storage<
          dynamical_system, symbol<"vd">,
          some::vdescriptor<typename dynamical_system_graphs::type> >,
      attached_storage<
          interaction, symbol<"vd">,
          some::vdescriptor<typename interaction_graphs::type> > >;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    decltype(auto) dynamical_system_graphs()
    {
      return Handle ::type ::dynamical_system_graphs ::at(*self());
    };
    decltype(auto) interaction_graphs()
    {
      return Handle ::type ::interaction_graphs ::at(*self());
    };
    decltype(auto) nids() { return Handle ::type ::nids ::at(*self()); };

    using default_interface<Handle>::self;

    template <match::handle<dynamical_system> Hds>
    decltype(auto) link(Hds ds)
    {
      auto &data = self()->data();
      auto &dsg0 = self()->dynamical_system_graphs()[0];
      auto &ig0 = self()->interaction_graphs()[0];
      ;

      auto inter = add<interaction>(data);
      auto dsgv = dsg0.add_vertex(ds);
      auto [dsg0_ed, ig0_vd] = dsg0.add_edge(dsgv, dsgv, inter, ig0);
      ds.property(symbol<"vd">{}) = dsgv;
      inter.property(symbol<"vd">{}) = ig0_vd;
      inter.property(symbol<"nds">{}) = 1;
      inter.template property<"ds1">() = ds;

      return inter;
    };

    template <match::handle<dynamical_system> Hds>
    decltype(auto) link(Hds ds1, Hds ds2)
    {
      auto &data = self()->data();
      auto &dsg0 = self()->dynamical_system_graphs()[0];
      auto &ig0 = self()->interaction_graphs()[0];
      ;

      auto inter = add<interaction>(data);
      auto dsgv1 = dsg0.add_vertex(ds1);
      auto dsgv2 = dsg0.add_vertex(ds2);
      auto [dsg0_ed, ig0_vd] = dsg0.add_edge(dsgv1, dsgv2, inter, ig0);

      ds1.property(symbol<"vd">{}) = dsgv1;
      ds2.property(symbol<"vd">{}) = dsgv2;

      inter.property(symbol<"nds">{}) = 2;
      inter.property(symbol<"vd">{}) = ig0_vd;
      inter.property(symbol<"ds1">{}) = ds1;
      inter.property(symbol<"ds2">{}) = ds2;

      return inter;
    };

    // strategy 2 : h_matrix is assembled for involved ds only
    auto make_index()
    {
      auto &data = self()->data();
      using info_t = std::decay_t<decltype(ground::get<info>(data))>;
      using env = typename info_t::env;
      using indice = typename env::indice;
      auto &dsg0 = self()->dynamical_system_graphs()[0];
      indice counter = 0;
      for (auto [dsi, dsiend] = dsg0.vertices(); dsi != dsiend; ++dsi) {
        indice &prop =
            handle(dsg0.bundle(*dsi), data).property(symbol<"index">{});
        auto &&involved =
            handle(dsg0.bundle(*dsi), data).property(symbol<"involved">{});
        auto [oei, oeiend] = dsg0.out_edges(*dsi);
        if (oeiend != oei)  // ds is involved in some interaction
        {
          involved = true;
          prop = counter++;
        }
        else {
          involved = false;
          prop = 0;
        }
      }
      self()->nids() = counter;
    };
  };
};
}  // namespace siconos
