#pragma once

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <typeinfo>

#include "siconos/collision/diskdisk_r.hpp"
#include "siconos/collision/diskfdisk_r.hpp"
#include "siconos/collision/disksegment_r.hpp"
#include "siconos/collision/point.hpp"
#include "siconos/collision/shape/disk.hpp"
#include "siconos/collision/shape/segment.hpp"
#include "siconos/collision/space_filter.hpp"
#include "siconos/siconos.hpp"
#include "siconos/storage/ground/ground.hpp"
#include "siconos/storage/pattern/pattern.hpp"
#include "siconos/storage/storage.hpp"

namespace py = pybind11;

namespace siconos::config::disks {
using disk = model::lagrangian_ds;
using nslaw = model::newton_impact_friction;
using diskdisk_r = collision::diskdisk_r;
using diskfdisk_r = collision::diskfdisk_r;
using disksegment_r = collision::disksegment_r;
using segment_shape = collision::shape::segment;
using disk_shape = collision::shape::disk;
using translated_disk_shape = collision::translated<disk_shape>;

using fc2d = simul::nonsmooth_problem<FrictionContactProblem>;
using osnspb = simul::one_step_nonsmooth_problem<fc2d>;
using solver_options = simul::solver_options;
using interaction =
    simul::interaction<nslaw, diskdisk_r, diskfdisk_r, disksegment_r>;
using osi = simul::one_step_integrator<disk, interaction>::moreau_jean;
using td = simul::time_discretization<>;
using topo = simul::topology<disk, interaction>;
using pointd = collision::point<disk>;
using pointl = collision::point<collision::shape::segment>;
using pointtds = collision::point<translated_disk_shape>;
using neighborhood = collision::neighborhood<pointd, pointl, pointtds>;
using space_filter = collision::space_filter<topo, neighborhood>;
using interaction_manager = simul::interaction_manager<space_filter>;
using simulation = simul::time_stepping<td, osi, osnspb, topo>;

using io = io::io<osi>;
using params = map<iparam<"dof", 3>, iparam<"ncgroups", 1>>;
}  // namespace siconos::config::disks

using namespace siconos;
namespace some = siconos::storage::some;

namespace pattern = siconos::storage::pattern;

namespace siconos::python::disks {

namespace config = siconos::config::disks;

static auto imake_storage()
{
  return storage::make<
      standard_environment<config::params>, config::simulation,
      config::interaction_manager, config::neighborhood, config::space_filter,
      config::io, config::disk, config::diskdisk_r, config::disksegment_r,
      config::diskfdisk_r, config::pointl, config::pointd, config::pointtds,
      config::interaction, config::segment_shape, config::disk_shape,
      storage::with_properties<
          storage::wrapped<config::disk, some::unbounded_collection>,
          storage::wrapped<config::diskdisk_r, some::unbounded_collection>,
          storage::wrapped<config::disksegment_r, some::unbounded_collection>,
          storage::wrapped<config::diskfdisk_r, some::unbounded_collection>,
          storage::wrapped<config::pointl, some::unbounded_collection>,
          storage::wrapped<config::pointd, some::unbounded_collection>,
          storage::wrapped<config::pointtds, some::unbounded_collection>,
          storage::wrapped<config::interaction, some::unbounded_collection>,
          storage::wrapped<config::segment_shape, some::unbounded_collection>,
          storage::wrapped<config::disk_shape, some::unbounded_collection>,
          storage::wrapped<config::translated_disk_shape,
                           some::unbounded_collection>,
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
          storage::bind<config::diskfdisk_r, "diskfdisk_r">,
          storage::bind<config::disksegment_r, "disksegment_r">,
          storage::bind<config::neighborhood, "neighborhood">,
          storage::bind<config::space_filter, "space_filter">,
          storage::bind<config::segment_shape, "segment_shape">,
          storage::bind<config::disk_shape, "disk_shape">,
          storage::bind<config::translated_disk_shape,
                        "translated_disk_shape">,
          storage::bind<config::interaction_manager, "interaction_manager">,
          storage::bind<config::interaction, "interaction">,
          storage::bind<config::osnspb, "osnspb">,
          storage::bind<config::solver_options, "solver_options">,
          storage::bind<config::osi, "osi">,
          storage::bind<config::td, "time_discretization">,
          storage::bind<config::topo, "topology">,
          storage::bind<config::simulation, "simulation">,
          storage::bind<config::osnspb, "osnspb">,
          storage::bind<config::fc2d, "fc2d">,
          storage::bind<config::io, "io">>>();
}

static auto idata = imake_storage();

using idata_t = std::decay_t<decltype(idata)>;

// static_assert(storage::pattern::match::diagonal_matrix<
//               decltype(storage::get<
//                        storage::pattern::attr_t<config::disk,
//                        "mass_matrix">>(
//                   idata_t{}, 0, storage::add<config::disk>(idata_t{})))>);

// just hide idata_t to pybind11
struct data_t {
  data_t() : _data(idata_t{}){};

  idata_t& operator()() { return _data; };

  idata_t _data;
};

}  // namespace siconos::python::disks

// namespace PYBIND11_NAMESPACE {
// namespace detail {
// template <>
// struct type_caster<siconos::python::disks::data_t> {};
// }}

namespace ground = siconos::storage::ground;
namespace match = siconos::storage::pattern::match;
template <typename H, typename T>
static decltype(auto) out_formatter(H h, T&& out_value)
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
static decltype(auto) in_formatter(H&& h, T&& in_value)
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

// static auto py_handles(py::module& mod)
// {
//   using disks_info_t = std::decay_t<decltype(ground::get<storage::info>(
//       siconos::python::disks::idata_t{}))>;

//   using disks_properties_t = typename disks_info_t::all_properties_t;

//   using disks_items_t = decltype(ground::transform(
//       typename disks_info_t::all_items_t{}, []<match::item I>(I) {
//         if constexpr (match::wrap<I>) {
//           return typename I::type{};
//         }
//         else {
//           return I{};
//         }
//       }));

//   // ground::type_trace<disks_items_t>();
//   auto named_disks_items = ground::filter(
//       disks_items_t{}, ground::is_a_model<[]<typename T>() {
//         return storage::has_property_from<T, storage::property::bind,
//                                           disks_properties_t>();
//       }>);

//   // ground::type_trace<std::decay_t<decltype(named_disks_items)>>();

//   auto disks_handles = ground::transform(
//       // only named items
//       named_disks_items, []<match::item I>(I item) {
//         // get handle
//         return storage::add<I>(siconos::python::disks::idata_t{});
//       });

//   // add corresponding py::class_
//   auto pyhandles = ground::transform(disks_handles, [&mod]<typename H>(
//                                                         H handle) {
//     using item_t = typename H::type;
//     using base_index_t = typename H::base_index_t;
//     auto base_index = py::class_<base_index_t>(
//         mod, fmt::format("index_{}",
//                          storage::bind_name<item_t, disks_properties_t>())
//                  .c_str());
//     return ground::make_tuple(
//         base_index,
//         py::class_<H>(mod, storage::bind_name<item_t,
//         disks_properties_t>(),
//                       base_index),
//         ground::type_c<H>);
//   });

//   return pyhandles;
// }
