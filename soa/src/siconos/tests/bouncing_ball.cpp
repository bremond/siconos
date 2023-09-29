
#include "siconos/siconos.hpp"
#include "siconos/utils/print.hpp"

namespace siconos::config {
using ball = model::lagrangian_ds;
using lcp = simul::nonsmooth_problem<LinearComplementarityProblem>;
using osnspb = simul::one_step_nonsmooth_problem<lcp>;
using nslaw = model::newton_impact;
using relation = model::lagrangian_tir<nslaw>;
using interaction = simul::interaction<relation>;
using osi = simul::one_step_integrator<ball, interaction>::moreau_jean;
using td = simul::time_discretization<>;
using topo = simul::topology<ball, interaction>;
using simulation = simul::time_stepping<td, osi, osnspb, topo>;

using params = map<iparam<"dof", 3>>;
}  // namespace siconos::config

template <typename E, typename C>
struct concat_env : E, C {
};

int main(int argc, char* argv[])
{
  using namespace siconos;
  auto data = storage::make_storage<
      standard_environment<config::params>, config::simulation, config::ball,
      config::relation, config::interaction,
      storage::with_properties<
          storage::diagonal<config::ball::mass_matrix>,
          storage::unbounded_diagonal<config::osi::mass_matrix_assembled>>>();

  // unsigned int nDof = 3;         // degrees of freedom for the ball
  double t0 = 0;               // initial computation time
  double tmax = 10;            // final computation time
  double h = 0.005;            // time step
  double position_init = 1.0;  // initial position for lowest bead.
  double velocity_init = 0.0;  // initial velocity for lowest bead.
  double theta = 0.5;          // theta for MoreauJeanOSI integrator
  double radius = 0.1;         // Ball radius
  double m = 1.;               // Ball mass
  double g = 9.81;             // Gravity

  print("====> Model loading ...\n");

  // --------------------------
  // -- The dynamical_system --
  // --------------------------
  auto ball = storage::add<config::ball>(data);

  ball.q() = {position_init, 0, 0};
  ball.velocity() = {velocity_init, 0, 0};
  ball.mass_matrix().diagonal() << m, m, 2. / 5. * m * radius * radius;

  // -- Set external forces (weight) --
  ball.fext() = {-m * g, 0., 0.};

  // ------------------
  // -- The relation --
  // ------------------

  // -- Lagrangian relation --
  auto relation = storage::add<config::relation>(data);
  //  the_relation.h_matrix() = {1.0, 0., 0.};

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

  // Interaction ball-floor
  auto interaction = simulation.topology().link(ball);
  interaction.h_matrix1() = {1., 0., 0.};

  interaction.relation() = relation;
  interaction.nonsmooth_law() = nslaw;

  // =========================== End of model definition
  // ===========================
  // ================================= Computation
  // =================================

  //  auto fd = io::open<ascii>("result.dat");

  // fix this for constant fext
  simulation.one_step_integrator().compute_iteration_matrix(
      simulation.current_step());

  auto out = fmt::output_file("result.dat");

  out.print("{:.15e} {:.15e} {:.15e} {:.15e} {:.15e}\n",
            simulation.current_step() * simulation.time_step(),
            config::ball::q::at(ball, simulation.current_step())(0),
            config::ball::velocity::at(ball, simulation.current_step())(0),
            0., 0.);

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
              config::ball::q::at(ball, simulation.current_step())(0),
              config::ball::velocity::at(ball, simulation.current_step())(0),
              p0, lambda);
  }
  //  io::close(fd);

}
