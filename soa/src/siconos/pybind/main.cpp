#include <pybind11/pybind11.h>

#include "siconos/siconos.hpp"

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
namespace py = pybind11;
namespace some = siconos::storage::some;
namespace config = siconos::config::disks;
using siconos::storage::pattern::wrap;

auto make_storage()
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

using data_t = std::decay_t<decltype(make_storage())>;

auto add_disk(data_t& data) { return storage::add<config::disk>(data); };
}  // namespace siconos::python::disks

PYBIND11_MODULE(nonos, m)
{
  auto disks = m.def_submodule("disks");
  disks.doc() = R"pbdoc(
        Nonos disks
        -----------

        .. currentmodule:: nonos

        .. autosummary::
           :toctree: _generate

           add
           subtract
    )pbdoc";

  disks.def("make_storage", &siconos::python::disks::make_storage, R"pbdoc(
        Create a new data object for 2D disks simulation
    )pbdoc");

  disks.def("add_disk", &siconos::python::disks::add_disk, R"pbdoc(
        Add disk
    )pbdoc");

  disks.attr("__version__") = "dev";
}
