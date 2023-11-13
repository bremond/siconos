#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include "siconos/siconos.hpp"

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

using params = map<iparam<"dof", 3>>;
}  // namespace siconos::config::disks

namespace siconos::python::disks {

using namespace siconos;
namespace some = siconos::storage::some;
namespace config = siconos::config::disks;
using siconos::storage::pattern::wrap;

auto imake_storage()
{
  return storage::make<
      standard_environment<config::params>, config::simulation,
      wrap<some::unbounded_collection, config::disk>, config::diskdisk_r,
      config::diskplan_r,
      wrap<some::unbounded_collection, config::interaction>,
      storage::with_properties<
          storage::attached<config::disk, storage::pattern::symbol<"shape">,
                            storage::some::item_ref<model::disk>>,
          storage::time_invariant<config::disk::fext>,
          storage::diagonal<config::disk, "mass_matrix">,
          storage::unbounded_diagonal<config::osi::mass_matrix_assembled>>>();
}

using idata_t = std::decay_t<decltype(imake_storage())>;

struct data_t {
  data_t() : _data(imake_storage()){};

  idata_t& operator()() { return _data; };

  idata_t _data;
};

data_t make_storage() { return data_t(); };

decltype(auto) add_disk(data_t& data)
{
  return storage::add<config::disk>(data());
};

using disk_t = std::decay_t<decltype(storage::add<config::disk>(idata_t{}))>;

}  // namespace siconos::python::disks

// namespace PYBIND11_NAMESPACE {
// namespace detail {
// template <>
// struct type_caster<siconos::python::disks::data_t> {};
// }}
PYBIND11_MODULE(_nonos, m)
{
  using disk_t = siconos::python::disks::disk_t;
  auto disks = m.def_submodule("disks");
  py::class_<siconos::python::disks::data_t>(disks, "data_t");
  py::class_<disk_t>(disks, "disk_t")
      .def("q", &disk_t::q, py::return_value_policy::reference)
      .def("velocity", &disk_t::velocity, py::return_value_policy::reference)
      .def("mass_matrix", &disk_t::mass_matrix,
           py::return_value_policy::reference);

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

  disks.def("add_disk", &siconos::python::disks::add_disk,
            py::return_value_policy::reference, R"pbdoc(
        Add disk
    )pbdoc");

  disks.attr("__version__") = "dev";
  m.attr("__version__") = "dev";
}
