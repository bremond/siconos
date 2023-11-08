#include <pybind11/pybind11.h>

#include "siconos/siconos.hpp"

namespace siconos::config {

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
}  // namespace siconos::config

using namespace siconos;
namespace py = pybind11;
namespace some = siconos::storage::some;
using siconos::storage::pattern::wrap;

auto new_2d_disks_data()
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

PYBIND11_MODULE(nonos, m)
{
  m.doc() = R"pbdoc(
        Nonos example
        -------------

        .. currentmodule:: nonos

        .. autosummary::
           :toctree: _generate

           add
           subtract
    )pbdoc";

  m.def("new_2d_disks_data", &new_2d_disks_data, R"pbdoc(
        Create a new data object for 2D disks simulation
    )pbdoc");

  m.attr("__version__") = "dev";
}
