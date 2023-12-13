#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <typeinfo>

#include "siconos/siconos.hpp"
#include "siconos/storage/ground/ground.hpp"
#include "siconos/storage/pattern/pattern.hpp"
#include "siconos/storage/storage.hpp"

namespace py = pybind11;

namespace siconos::config::disks {
using disk = model::lagrangian_ds;
using lcp = simul::nonsmooth_problem<LinearComplementarityProblem>;
using osnspb = simul::one_step_nonsmooth_problem<lcp>;
using nslaw = model::newton_impact;
using diskdisk_r = model::diskdisk_r;
using diskplan_r = model::diskplan_r;
using interaction = simul::interaction<nslaw, diskdisk_r, diskplan_r>;
using osi = simul::one_step_integrator<disk, interaction>::moreau_jean;
using td = simul::time_discretization<>;
using topo = simul::topology<disk, interaction>;
using simulation = simul::time_stepping<td, osi, osnspb, topo>;
using disk_shape = model::disk_shape;
using params = map<iparam<"dof", 3>>;
}  // namespace siconos::config::disks

using namespace siconos;
namespace some = siconos::storage::some;

namespace pattern = siconos::storage::pattern;

namespace siconos::python::disks {

namespace config = siconos::config::disks;

auto imake_storage()
{
  return storage::make<
      standard_environment<config::params>, config::simulation,
      pattern::wrap<some::unbounded_collection, config::disk>,
      pattern::wrap<some::unbounded_collection, config::diskdisk_r>,
      pattern::wrap<some::unbounded_collection, config::diskplan_r>,
      pattern::wrap<some::unbounded_collection, config::interaction>,
      pattern::wrap<some::unbounded_collection, config::disk_shape>,
      storage::with_properties<
          storage::attached<config::disk, storage::pattern::symbol<"shape">,
                            storage::some::item_ref<config::disk_shape>>,
          storage::time_invariant<
              storage::pattern::attr_t<config::disk, "fext">>,
          storage::diagonal<
              storage::pattern::attr_t<config::disk, "mass_matrix">>,
          storage::unbounded_diagonal<
              storage::pattern::attr_t<config::osi, "mass_matrix_assembled">>,
          storage::bind<config::disk, "disk">,
          storage::bind<config::nslaw, "nslaw">,
          storage::bind<config::diskdisk_r, "diskdisk_r">,
          storage::bind<config::diskplan_r, "diskplan_r">,
          storage::bind<config::disk_shape, "disk_shape">,
          storage::bind<config::interaction, "interaction">,
          storage::bind<config::osi, "osi">,
          storage::bind<config::td, "time_discretization">,
          storage::bind<config::topo, "topology">,
          storage::bind<config::simulation, "simulation">,
          storage::bind<config::osnspb, "osnspb">,
          storage::bind<config::lcp, "lcp">>>();
}

using idata_t = std::decay_t<decltype(imake_storage())>;

static_assert(storage::pattern::match::diagonal_matrix<
              decltype(storage::get<
                       storage::pattern::attr_t<config::disk, "mass_matrix">>(
                  idata_t{}, 0, storage::add<config::disk>(idata_t{})))>);

// just hide idata_t to pybind11
struct data_t {
  data_t() : _data(imake_storage()){};

  idata_t& operator()() { return _data; };

  idata_t _data;
};

data_t make_storage() { return data_t(); };

}  // namespace siconos::python::disks

// namespace PYBIND11_NAMESPACE {
// namespace detail {
// template <>
// struct type_caster<siconos::python::disks::data_t> {};
// }}

namespace ground = siconos::storage::ground;
namespace match = siconos::storage::pattern::match;
template <typename H, typename T>
decltype(auto) out_formatter(H h, T&& out_value)
{
  using out_t = std::decay_t<T>;

  if constexpr (!match::diagonal_matrix<out_t>) {
    if constexpr (match::matrix<out_t>) {
      return algebra::matrix_ref<out_t>(static_cast<T&&>(out_value));
    }
    else if constexpr (match::index<out_t>) {
      return storage::handle(h.data(), static_cast<T&&>(out_value));
    }
    else {
      return static_cast<T&&>(out_value);
    }
  }
  else {
    return out_value.diagonal();
  }
}

template <typename H, typename T>
decltype(auto) in_formatter(H&& h, T&& in_value)
{
  using in_t = std::decay_t<T>;

  if constexpr (!match::diagonal_matrix<in_t>) {
    return static_cast<T&&>(in_value);
  }
  else {
    return in_value.diagonal();
  }
}

using namespace boost::hana::literals;
PYBIND11_MODULE(_nonos, m)
{
  // sub module disks
  auto disks = m.def_submodule("disks");

  auto data_class =
      py::class_<siconos::python::disks::data_t>(disks, "data_t");

  using disks_info_t = std::decay_t<decltype(ground::get<storage::info>(
      siconos::python::disks::idata_t{}))>;

  using disks_properties_t = typename disks_info_t::all_properties_t;

  using disks_items_t =
    decltype(
      ground::transform(typename disks_info_t::all_items_t{},
                        []<match::item I>(I) {
                          if constexpr (match::wrap<I>) {
                            return typename I::type{};
                          }
                          else {
                            return I{};
                          }
                        }));

      //  ground::type_trace<disks_items_t>();
      auto named_disks_items = ground::filter(
          disks_items_t{}, ground::is_a_model<[]<typename T>() {
            return storage::has_property_from<T, storage::property::bind,
                                              disks_properties_t>();
          }>);

  auto disks_handles = ground::transform(
      // only named items
      named_disks_items, []<match::item I>(I item) {
        // get handle
        return storage::add<I>(siconos::python::disks::idata_t{});
      });

  auto pyhandles =
      ground::transform(disks_handles, [&disks]<typename H>(H handle) {
        using item_t = typename H::type;
        return ground::make_tuple(
            py::class_<H>(disks,
                          storage::bind_name<item_t, disks_properties_t>()),
            ground::type_c<H>);
      });

  // attached storage
  ground::for_each(pyhandles, [](auto pyhandle) {
    using handle_t = typename decltype(+pyhandle[1_c])::type;
    using item_t = typename handle_t::type;

    ground::fold_left(
        decltype(storage::attached_storages(
            handle_t{}, handle_t{}.data())){},  // all attached storages
        std::ref(pyhandle[0_c]),                // initial state
        []<match::attached_storage<item_t> S>(py::class_<handle_t> dc, S s) {
          constexpr auto astor_name = storage::attached_storage_name(s);
          using target_type = std::decay_t<decltype(out_formatter(
              handle_t{}, storage::prop<astor_name.str>(handle_t{})))>;

          return dc
              .def(
                  fmt::format("{}", astor_name.str.value).c_str(),
                  [&astor_name](handle_t& h) -> target_type {
                    return out_formatter(h, storage::prop<astor_name.str>(h));
                  },
                  py::return_value_policy::reference)
              .def(fmt::format("set_{}", astor_name.str.value).c_str(),
                   [&astor_name](handle_t& h, target_type val) {
                     in_formatter(h, storage::prop<astor_name.str>(h)) =
                         in_formatter(h, val);
                   });
        });
  });

  // attributes
  ground::for_each(pyhandles, [](auto pyhandle) {
    using handle_t = typename decltype(+pyhandle[1_c])::type;
    //    using item_t = typename handle_t::type;
    ground::fold_left(
        pattern::attributes(typename handle_t::type{}),  // all attributes
        std::ref(pyhandle[0_c]),                         // initial state
        []<match::attribute A>(py::class_<handle_t> dc, A a) {
          using attr_value_t = std::decay_t<decltype(out_formatter(
              handle_t{},
              storage::get<A>(handle_t{}.data(), 0, handle_t{})))>;
          return dc
              .def(
                  fmt::format("{}", pattern::attribute_name(a)).c_str(),
                  [](handle_t& h) -> attr_value_t {
                    return out_formatter(h, storage::get<A>(h.data(), h));
                  },
                  py::return_value_policy::reference)
              .def(fmt::format("set_{}", pattern::attribute_name(a)).c_str(),
                   [](handle_t& h, attr_value_t v) {
                     in_formatter(h, storage::get<A>(h.data(), h)) =
                         in_formatter(h, v);
                   });
        });
  });

  // methods
  ground::for_each(pyhandles, [](auto pyhandle) {
    using handle_t = typename decltype(+pyhandle[1_c])::type;
    //    using item_t = typename handle_t::type;
    ground::fold_left(storage::methods(pyhandle[1_c]),  // all methods
                      std::ref(pyhandle[0_c]),          // initial state
                      []<typename M>(py::class_<handle_t> dc, M m) {
                        return dc.def(pattern::method_name(m),
                                      pattern::method_def(m),
                                      py::return_value_policy::reference);
                      });
  });

  disks.doc() = R"pbdoc(
        Nonos m
        -----------

        .. currentmodule:: nonos

        .. autosummary::
           :toctree: _generate

           add
           subtract
    )pbdoc";

  disks.def("make_storage", &siconos::python::disks::make_storage,
            py::return_value_policy::reference,
            R"pbdoc(
        Create a new data object for 2D disks simulation
    )pbdoc");

  ground::for_each(pyhandles, [&disks](auto pyhandle) {
    using handle_t = typename decltype(+pyhandle[1_c])::type;
    using item_t = typename handle_t::type;
    auto item_name = storage::bind_name<item_t, disks_properties_t>();

    disks.def(
        fmt::format("add_{}", item_name).c_str(),
        [](siconos::python::disks::data_t& data) {
          return siconos::storage::add<item_t>(data());
        },
        py::return_value_policy::reference, R"pbdoc(
        Add
    )pbdoc");
  });

  disks.attr("__version__") = "dev";
  m.attr("__version__") = "dev";
}
