
#include "siconos/siconos.hpp"
#include "siconos/utils/print.hpp"

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

int main(int argc, char* argv[])
{
  using namespace siconos;
  auto data = storage::make<
      standard_environment<config::params>, config::simulation, config::disk,
      config::diskdisk_r, config::diskplan_r, config::interaction,
      storage::with_properties<
          storage::attached<config::disk, storage::pattern::symbol<"shape">,
                            storage::some::item_ref<model::disk>>,
          storage::time_invariant<storage::attr_of<config::disk, "fext">>,
          storage::diagonal<config::disk, "mass_matrix">,
          storage::unbounded_diagonal<config::osi::mass_matrix_assembled>>>();

  // unsigned int nDof = 3;         // degrees of freedom for the disk
  double t0 = 0;               // initial computation time
  double tmax = 10;            // final computation time
  double h = 0.005;            // time step
  double position_init = 1.0;  // initial position for lowest bead.
  double velocity_init = 0.0;  // initial velocity for lowest bead.
  double theta = 0.5;          // theta for MoreauJeanOSI integrator
  double radius = 0.1;         // Disk radius
  double m = 1.;               // Disk mass
  double g = 9.81;             // Gravity

  print("====> Model loading ...\n");

  // --------------------------
  // -- The dynamical_system --
  // --------------------------
  auto d1 = storage::add<config::disk>(data);
  auto d2 = storage::add<config::disk>(data);

  d1.q() = {position_init, 0, 0};
  d1.velocity() = {velocity_init, 0, 0};
  d1.mass_matrix().diagonal() << m, m, 2. / 5. * m * radius * radius;

  d2.q() = {2 * position_init, 0.1, 0};
  d2.velocity() = {velocity_init, 0, 0};
  d2.mass_matrix().diagonal() << m, m, 2. / 5. * m * radius * radius;

  storage::handle(data, prop<"shape">(d1)).radius() = 1;
  storage::handle(data, prop<"shape">(d2)).radius() = 2;

  // -- Set external forces (weight) --
  d1.fext() = {-m * g, 0., 0.};
  d2.fext() = {-m * g, 0., 0.};

  // ------------------
  // -- The relation --
  // ------------------

  // -- Lagrangian relation --
  auto dd_r = storage::add<config::diskdisk_r>(data);
  auto d1p_r = storage::add<config::diskplan_r>(data);
  auto d2p_r = storage::add<config::diskplan_r>(data);

  // -- nslaw --
  double e = 0.9;
  auto nslaw = storage::add<config::nslaw>(data);
  nslaw.e() = e;

  auto lcp = storage::add<config::lcp>(data);
  lcp.create();

  // ------------------
  // --- Simulation ---
  // ------------------
  auto simulation = storage::add<config::simulation>(data);

  simulation.one_step_integrator().theta() = theta;
  simulation.one_step_integrator().constraint_activation_threshold() = 0.;
  simulation.time_discretization().t0() = t0;
  simulation.time_discretization().h() = h;

  simulation.time_discretization().tmax() = tmax;
  // -- set the formulation for the one step nonsmooth problem --
  auto osnspb = simulation.one_step_nonsmooth_problem();
  osnspb.problem() = lcp;

  // -- set the options --
  auto so = storage::add<simul::solver_options>(data);
  so.create();
  osnspb.options() = so;

  // Interaction disk-disk
  auto inter_dd = simulation.topology().link(d1, d2);

  // Interaction disk-plan
  auto inter_d1p = simulation.topology().link(d1);
  auto inter_d2p = simulation.topology().link(d2);

  inter_dd.relation() = dd_r;
  inter_dd.nonsmooth_law() = nslaw;

  inter_d1p.relation() = d1p_r;
  inter_d1p.nonsmooth_law() = nslaw;

  inter_d2p.relation() = d2p_r;
  inter_d2p.nonsmooth_law() = nslaw;

  // =========================== End of model definition
  // ===========================
  // ================================= Computation
  // =================================

  //  auto fd = io::open<ascii>("result.dat");

  // fix this for constant fext
  simulation.initialize();

  auto out = fmt::output_file("result.dat");

  out.print("{:.15e} {:.15e} {:.15e} {:.15e} {:.15e}\n",
            simulation.current_step() * simulation.time_step(),
            storage::attr<"q">(d1, simulation.current_step())(0),
            storage::attr<"velocity">(d1, simulation.current_step())(0), 0.,
            0.);

  while (simulation.has_next_event()) {
    auto ninvds = simulation.compute_one_step();

    double p0, lambda;
    if (ninvds > 0) {
      p0 = get_vector(simulation.one_step_integrator().p0_vector_assembled(),
                      0)(0);
      lambda = get_vector(
          simulation.one_step_integrator().lambda_vector_assembled(), 0)(0);
    }
    else {
      p0 = 0;
      lambda = 0;
    }

    out.print("{:.15e} {:.15e} {:.15e} {:.15e} {:.15e}\n",
              simulation.current_step() * simulation.time_step(),
              storage::attr<"q">(d1, simulation.current_step())(0),
              storage::attr<"velocity">(d1, simulation.current_step())(0), p0,
              lambda);
  }
  //  io::close(fd);
}
