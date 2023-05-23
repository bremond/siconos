#include <fmt/core.h>
#include <fmt/ranges.h>

#include "siconos/siconos.hpp"

using namespace siconos;
using env = standard_environment;

using fmt::print;
int main(int argc, char* argv[])
{
  using ball = lagrangian_ds;
  using lcp = numerics::nonsmooth_problem<LinearComplementarityProblem>;
  using osnspb = numerics::one_step_nonsmooth_problem<lcp>;
  using relation = lagrangian_r;
  using nslaw = nonsmooth_law::newton_impact;
  using interaction = interaction<nslaw, relation, 1>;
  using osi = one_step_integrator<ball, interaction>::moreau_jean;
  using td = time_discretization<>;
  using topo = topology<ball, interaction>;
  using simulation = time_stepping<td, osi, osnspb, topo>;
  using siconos::get;

  auto data = make_storage<
      standard_environment, simulation, ball, relation, interaction,
      with_properties<diagonal<ball::mass_matrix>,
                      unbounded_diagonal<osi::mass_matrix_assembled>>>();

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
  auto the_ball = add<ball>(data);

  the_ball.q() = {position_init, 0, 0};
  the_ball.velocity() = {velocity_init, 0, 0};
  the_ball.mass_matrix().diagonal() << m, m, 2. / 5. * m * radius * radius;

  // -- Set external forces (weight) --
  the_ball.fext() = {-m * g, 0., 0.};

  // ------------------
  // -- The relation --
  // ------------------

  // -- Lagrangian relation --
  auto the_relation = add<relation>(data);
  //  the_relation.h_matrix() = {1.0, 0., 0.};

  // -- nslaw --
  double e = 0.9;
  auto the_nslaw = add<nslaw>(data);
  the_nslaw.e() = e;

  auto the_lcp = add<lcp>(data);
  the_lcp.create();

  // ------------------
  // --- Simulation ---
  // ------------------
  auto the_simulation = add<simulation>(data);

  the_simulation.one_step_integrator().constraint_activation_threshold() = 0.;
  the_simulation.time_discretization().t0() = t0;
  the_simulation.time_discretization().h() = h;

  the_simulation.time_discretization().tmax() = tmax;
  // -- set the formulation for the one step nonsmooth problem --
  auto the_osnspb = the_simulation.one_step_nonsmooth_problem();
  the_osnspb.problem() = the_lcp;

  // -- set the options --
  auto so = add<numerics::solver_options>(data);
  so.create();
  the_osnspb.options() = so;

  // Interaction ball-floor
  auto the_interaction = the_simulation.topology().link(the_ball);
  the_interaction.h_matrix() = {1., 0., 0.};
  the_interaction.relation() = the_relation;
  the_interaction.nonsmooth_law() = the_nslaw;

  // =========================== End of model definition
  // ===========================
  // ================================= Computation
  // =================================

  //  auto fd = io::open<ascii>("result.dat");

  // fix this for constant fext
  the_simulation.one_step_integrator().compute_iteration_matrix(
      the_simulation.current_step());

  auto out = fmt::output_file("result.dat");
  while (the_simulation.has_next_event()) {
    the_simulation.compute_one_step();

    double p0, lambda;
    if (the_simulation.topology().ninvds() > 0) {
      p0 = get_vector(
          the_simulation.one_step_integrator().p0_vector_assembled(), 0)(0);
      lambda = get_vector(
          the_simulation.one_step_integrator().lambda_vector_assembled(),
          0)(0);
    }
    else {
      p0 = 0;
      lambda = 0;
    }

    out.print("{} {} {} {} {}\n",
              the_simulation.current_step() * the_simulation.time_step(),
              ball::q::at(the_ball, the_simulation.current_step())(0),
              ball::velocity::at(the_ball, the_simulation.current_step())(0),
              p0, lambda);
  }
  //  io::close(fd);
}
