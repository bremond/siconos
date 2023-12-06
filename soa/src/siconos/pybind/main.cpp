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
using siconos::storage::pattern::with_name;
using disk = with_name<"disk", model::lagrangian_ds>;
using lcp =
    with_name<"lcp", simul::nonsmooth_problem<LinearComplementarityProblem>>;
using osnspb = with_name<"osnspb", simul::one_step_nonsmooth_problem<lcp>>;
using nslaw = with_name<"nslaw", model::newton_impact>;
using diskdisk_r = with_name<"diskdisk_r", model::diskdisk_r>;
using diskplan_r = with_name<"diskplan_r", model::diskplan_r>;
using interaction =
    with_name<"interaction",
              simul::interaction<nslaw, diskdisk_r, diskplan_r>>;
using osi =
    with_name<"osi",
              simul::one_step_integrator<disk, interaction>::moreau_jean>;
using td = with_name<"time_discretization", simul::time_discretization<>>;
using topo = with_name<"topology", simul::topology<disk, interaction>>;
using simulation =
    with_name<"simulation", simul::time_stepping<td, osi, osnspb, topo>>;

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
      config::diskdisk_r, config::diskplan_r,
      pattern::wrap<some::unbounded_collection, config::interaction>,
      storage::with_properties<
          storage::attached<config::disk, storage::pattern::symbol<"shape">,
                            storage::some::item_ref<model::disk_shape>>,
          storage::time_invariant<
              storage::pattern::attr_t<config::disk, "fext">>,
          storage::diagonal<
              storage::pattern::attr_t<config::disk, "mass_matrix">>,
          storage::unbounded_diagonal<storage::pattern::attr_t<
              config::osi, "mass_matrix_assembled">>>>();
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

using disk_t = std::decay_t<decltype(storage::add<config::disk>(idata_t{}))>;

using nslaw_t =
    std::decay_t<decltype(storage::add<config::nslaw>(idata_t{}))>;

}  // namespace siconos::python::disks

// namespace PYBIND11_NAMESPACE {
// namespace detail {
// template <>
// struct type_caster<siconos::python::disks::data_t> {};
// }}

using namespace boost::hana::literals;
PYBIND11_MODULE(_nonos, m)
{
  namespace ground = siconos::storage::ground;
  namespace match = siconos::storage::pattern::match;

  // sub module disks
  auto disks = m.def_submodule("disks");

  auto data_class =
      py::class_<siconos::python::disks::data_t>(disks, "data_t");
  //  auto disk_class = py::class_<disk_t>(disks, "disk_t");
  //  auto nslaw_class = py::class_<nslaw_t>(disks, "nslaw_t");

  using disks_info_t = std::decay_t<decltype(ground::get<storage::info>(
      siconos::python::disks::idata_t{}))>;

  using disks_items_t = typename disks_info_t::all_items_t;

  auto disks_handles = ground::transform(
      // only named items
      ground::filter(disks_items_t{},
                     ground::derive_from<pattern::any_symbol>),
      []<match::item I>(I item) {
        // get handle
        return storage::add<I>(siconos::python::disks::idata_t{});
      });

  auto pyhandles =
      ground::transform(disks_handles, [&disks]<typename H>(H handle) {
        using item_t = typename H::type;
        return ground::make_tuple(
            py::class_<H>(disks, pattern::item_name(item_t{})),
            ground::type_c<H>);
      });

  ground::for_each(pyhandles, [](auto pyhandle) {
    using handle_t = typename decltype(+pyhandle[1_c])::type;
    using item_t = typename handle_t::type;
    ground::fold_left(
        ground::concat(
            decltype(storage::attached_storages(handle_t{},
                                                handle_t{}.data())){},
            ground::concat(pattern::attributes(typename handle_t::type{}),
                           storage::methods(pyhandle[1_c]))),
        std::ref(pyhandle[0_c]),
        ground::overload(
            []<match::attached_storage<item_t> S>(py::class_<handle_t> dc,
                                                  S s) {
              constexpr auto astor_name = storage::attached_storage_name(s);
              using target_type =
                  std::decay_t<decltype(storage::prop<astor_name.str>(
                      handle_t{}))>;

              return dc.def(
                  fmt::format("{}", astor_name.str.value).c_str(),
                  [&astor_name](handle_t& h) {
                    if constexpr (match::npy_format<target_type>) {
                      return storage::prop<astor_name.str>(h);
                    }
                    else {
                      // pseudo alert, fix with an exception ?
                      return std::move(
                          fmt::format(
                              "binding unimplemented for {}, type {} need "
                              "to be converted",
                              astor_name.str.value,
                              typeid(target_type).name())
                              .c_str());
                    }
                  },
                  py::return_value_policy::reference);
            },
            []<match::attribute A>(py::class_<handle_t> dc, A a) {
              using attr_value_t = std::decay_t<decltype(storage::get<A>(
                  handle_t{}.data(), 0, handle_t{}))>;
              if constexpr (!pattern::match::diagonal_matrix<attr_value_t>) {
                return dc
                    .def(
                        fmt::format("{}", pattern::attribute_name(a)).c_str(),
                        [](handle_t& h) {
                          if constexpr (match::index<std::decay_t<
                                            decltype(storage::get<A>(h.data(),
                                                                     h))>>) {
                            // return type is an item_ref: build a handle
                            return storage::handle(
                                h.data(), storage::get<A>(h.data(), h));
                          }
                          else {
                            // anything else
                            return storage::get<A>(h.data(), h);
                          }
                        },
                        py::return_value_policy::reference)
                    .def(fmt::format("set_{}", pattern::attribute_name(a))
                             .c_str(),
                         [](handle_t& h, attr_value_t v) {
                           storage::get<A>(h.data(), h) = v;
                         });
              }
              else {
                return dc;
              }
            },
            []<typename M>(py::class_<handle_t> dc, M m) {
              return dc.def(pattern::method_name(m), pattern::method_def(m));
            }));
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
            py::return_value_policy::reference, R"pbdoc(
        Create a new data object for 2D disks simulation
    )pbdoc");

  ground::for_each(pyhandles, [&disks](auto pyhandle) {
    using handle_t = typename decltype(+pyhandle[1_c])::type;
    using item_t = typename handle_t::type;
    auto item_name = pattern::item_name(item_t{});

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
