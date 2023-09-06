
#include "siconos/siconos.hpp"
#include "siconos/utils/print.hpp"

namespace siconos::data {
using ball = model::lagrangian_ds;
using lcp = numerics::nonsmooth_problem<LinearComplementarityProblem>;
using osnspb = numerics::one_step_nonsmooth_problem<lcp>;
using relation = model::lagrangian_tir;
using nslaw = model::nsl::newton_impact;
using interaction = simul::interaction<nslaw, relation, 1>;
using osi = simul::one_step_integrator<ball, interaction>::moreau_jean;
using td = simul::time_discretization<>;
using topo = simul::topology<ball, interaction>;
using simulation = simul::time_stepping<td, osi, osnspb, topo>;
}  // namespace siconos::data

int main(int argc, char* argv[])
{
  using namespace siconos;
  auto data = storage::make_storage<
      standard_environment, data::simulation,
      wrap<some::unbounded_collection, data::ball>,
    wrap<some::unbounded_collection, data::relation>,
    wrap<some::unbounded_collection, data::interaction>,
    storage::with_properties<
      storage::diagonal<data::ball::mass_matrix>,
      storage::unbounded_diagonal<data::osi::mass_matrix_assembled>>>();

  // unsigned int nDof = 3;         // degrees of freedom for the ball
  double t0 = 0;               // initial computation time
  double tmax = 1;             // final computation time
  double h = 0.005;            // time step
  double position_init = 1.0;  // initial position for lowest bead.
  double velocity_init = 0.0;  // initial velocity for lowest bead.
  double theta = 0.5;          // theta for MoreauJeanOSI integrator
  double radius = 0.1;         // Ball radius
  double m = 1.;               // Ball mass
  double g = 9.81;             // Gravity

  print("====> Model loading ...\n");

  // ---------------------------
  // -- The dynamical_systems --
  // ---------------------------
  auto ball1 = storage::add<data::ball>(data);
  auto ball2 = storage::add<data::ball>(data);

  ball1.q() = {position_init, 0, 0};
  ball2.q() = {position_init * 2, 0, 0};

  fmt::print("{} - {}\n", ball1.q(), ball2.q());

  ball1.velocity() = {velocity_init, 0, 0};
  ball2.velocity() = {velocity_init, 0, 0};
  ball1.mass_matrix().diagonal() << m, m, 2. / 5. * m * radius * radius;
  ball2.mass_matrix().diagonal() << m, m, 2. / 5. * m * radius * radius;

  // -- Set external forces (weight) --
  ball1.fext() = {-m * g, 0., 0.};
  ball2.fext() = {-m * g, 0., 0.};

  // ------------------
  // -- The relation --
  // ------------------

  // -- Lagrangian relation --
  auto relation_f = storage::add<data::relation>(data);
  auto relation_b = storage::add<data::relation>(data);
  relation_f.b() = - radius;
  relation_b.b() = -2 * radius;

  auto vbs = storage::attr_values<typename data::relation::b>(data, 0);
//  assert(vbs[0] == 0);
//  assert(vbs[1] == -0.2);
  //  the_relation.h_matrix() = {1.0, 0., 0.};

  // -- nslaw --
  double e = 0.9;
  auto nslaw = storage::add<data::nslaw>(data);
  nslaw.e() = e;

  auto lcp = storage::add<data::lcp>(data);
  lcp.create();

  // ------------------
  // --- Simulation ---
  // ------------------
  auto simulation = storage::add<data::simulation>(data);

  simulation.one_step_integrator().theta() = theta;
  simulation.one_step_integrator().constraint_activation_threshold() = 0.;
  simulation.time_discretization().t0() = t0;
  simulation.time_discretization().h() = h;

  simulation.time_discretization().tmax() = tmax;
  // -- set the formulation for the one step nonsmooth problem --
  auto osnspb = simulation.one_step_nonsmooth_problem();
  osnspb.problem() = lcp;

  // -- set the options --
  auto so = storage::add<numerics::solver_options>(data);
  so.create();
  osnspb.options() = so;

  // Interaction ball-floor
  auto interaction1 = simulation.topology().link(ball1);
  auto interaction2 = simulation.topology().link(ball1, ball2);
  interaction1.h_matrix1() = {1., 0., 0.};
  interaction1.h_matrix2() = {1., 0., 0.};
  interaction2.h_matrix1() = {-1., 0., 0.};
  interaction2.h_matrix2() = {1., 0., 0.};

  interaction1.relation() = relation_f;
  interaction2.relation() = relation_b;

  interaction1.nonsmooth_law() = nslaw;
  interaction2.nonsmooth_law() = nslaw;

  // =========================== End of model definition
  // ===========================
  // ================================= Computation
  // =================================

  //  auto fd = io::open<ascii>("result.dat");

  // fix this for constant fext
  simulation.one_step_integrator().compute_iteration_matrix(
      simulation.current_step());

  auto out = fmt::output_file("result-many.dat");

  out.print("{:.15e} {:.15e} {:.15e} {:.15e} {:.15e} {:.15e} {:.15e}\n",
            simulation.current_step() * simulation.time_step(),
            data::ball::q::at(ball1, simulation.current_step())(0),
            data::ball::q::at(ball2, simulation.current_step())(0),
            data::ball::velocity::at(ball1, simulation.current_step())(0),
            data::ball::velocity::at(ball2, simulation.current_step())(0), 0.,
            0.);

  while (simulation.has_next_event()) {
    auto ninvds = simulation.compute_one_step();

    double p01, p02, lambda1, lambda2;
    if (ninvds > 1) {
      p01 = get_vector(simulation.one_step_integrator().p0_vector_assembled(),
                       0)(0);
      p02 = get_vector(simulation.one_step_integrator().p0_vector_assembled(),
                       0)(1);

      lambda1 = get_vector(
          simulation.one_step_integrator().lambda_vector_assembled(), 0)(0);

      lambda2 = get_vector(
          simulation.one_step_integrator().lambda_vector_assembled(), 1)(0);
    }
    else if (ninvds == 1) {
      p01 = get_vector(simulation.one_step_integrator().p0_vector_assembled(),
                       0)(0);
      p02 = 0;

      lambda1 = get_vector(
          simulation.one_step_integrator().lambda_vector_assembled(), 0)(0);

      lambda2 = 0;
    }
    else {
      p01 = 0;
      p02 = 0;
      lambda1 = 0;
      lambda2 = 0;
    }

    out.print("{:.15e} {:.15e} {:.15e} {:.15e} {:.15e} {:.15e} {:.15e}\n",
              simulation.current_step() * simulation.time_step(),
              data::ball::q::at(ball1, simulation.current_step())(0),
              data::ball::q::at(ball2, simulation.current_step())(0),
              data::ball::velocity::at(ball1, simulation.current_step())(0),
              data::ball::velocity::at(ball2, simulation.current_step())(0),
              p01, p02, lambda1, lambda2);
  }
  //  io::close(fd);
}
