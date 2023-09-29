
#include <Friction_cst.h>
#include <stdlib.h>

#include <chrono>
#include <numeric>

#include "siconos/model/nslaws.hpp"
#include "siconos/siconos.hpp"
#include "siconos/utils/print.hpp"

namespace siconos::config {
using ball = model::lagrangian_ds;
using lcp = simul::nonsmooth_problem<LinearComplementarityProblem>;
using fc2d = simul::nonsmooth_problem<FrictionContactProblem>;
using osnspb = simul::one_step_nonsmooth_problem<fc2d>;
using nslaw = model::newton_impact_friction;
using relation = model::lagrangian_tir<nslaw>;
using interaction = simul::interaction<relation>;
using osi = simul::one_step_integrator<ball, interaction>::moreau_jean;
using td = simul::time_discretization<>;
using topo = simul::topology<ball, interaction>;
using simulation = simul::time_stepping<td, osi, osnspb, topo>;

using params = map<iparam<"dof", 3>>;

}  // namespace siconos::config

int main(int argc, char* argv[])
{
  unsigned int nballs = atoi(argv[1]);
  using namespace siconos;
  namespace some = siconos::storage::some;
  using siconos::storage::pattern::wrap;

  auto data = storage::make_storage<
      standard_environment<config::params>, config::simulation,
      wrap<some::unbounded_collection, config::ball>,
      wrap<some::unbounded_collection, config::relation>,
      wrap<some::unbounded_collection, config::interaction>,
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

  // ---------------------------
  // -- The dynamical_systems --
  // ---------------------------
  for (unsigned int i = 0; i < nballs; ++i) {
    auto ball = storage::add<config::ball>(data);
    ball.q() = {position_init * (i + 1), 0, 0};
    ball.velocity() = {velocity_init, 0, 0};
    ball.mass_matrix().diagonal() << m, m, 2. / 5. * m * radius * radius;
    ball.fext() = {-m * g, 0., 0.};
  }

  // for (auto ball : storage::handles<config::ball>(data, 0))
  //{
  //   print("ball:{} , ball.q()={}\n", ball.get(), ball.q()[0]);
  // }
  //  ------------------
  //  -- The relation --
  //  ------------------

  // -- Lagrangian relation --
  auto relation_f = storage::add<config::relation>(data);
  auto relation_b = storage::add<config::relation>(data);
  relation_f.b() = -radius;
  relation_b.b() = -2 * radius;

  //  assert(vbs[0] == 0);
  //  assert(vbs[1] == -0.2);
  //  the_relation.h_matrix() = {1.0, 0., 0.};

  // -- nslaw --
  double e = 0.9;
  auto nslaw = storage::add<config::nslaw>(data);
  nslaw.e() = e;
  nslaw.mu() = 0.;

  //  auto lcp = storage::add<config::lcp>(data);
  //  lcp.create();
  auto fc2d = storage::add<config::fc2d>(data);
  fc2d.create();
  fc2d.instance()->dimension = 2;
  fc2d.instance()->mu = 0;

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
  osnspb.problem() = fc2d;

  // -- set the options --
  auto so = storage::add<simul::solver_options>(data);
  so.create(SICONOS_FRICTION_2D_NSGS);
  so.instance()->iparam[SICONOS_IPARAM_MAX_ITER] = 1;
  osnspb.options() = so;

  auto balls = storage::handles<config::ball>(data, 0);

  auto first_ball = (balls | views::take(1)).front();
  //    views::transform([&simulation, &radius, &relation_f, &nslaw](auto
  //    first_ball)
  //    {
  auto interaction = simulation.topology().link(first_ball);
  interaction.h_matrix1() << 1., 0., 0., 0., 1., -radius;
  interaction.h_matrix2() << 1., 0., 0., 0., 1., -radius;
  interaction.relation() = relation_f;
  interaction.nonsmooth_law() = nslaw;
  //    });

  for (auto [ball1, ball2] : views::zip(balls, balls | views::drop(1))) {
    // print("new interaction ball<->ball : {} {}\n", ball1.get(),
    // ball2.get());
    auto interaction = simulation.topology().link(ball1, ball2);
    interaction.h_matrix1() << -1., 0., 0., 0., 1., -radius;
    interaction.h_matrix2() << 1., 0., 0., 0., 1., -radius;
    interaction.relation() = relation_b;
    interaction.nonsmooth_law() = nslaw;
  };

  // =========================== End of model definition
  // ===========================
  // ================================= Computation
  // =================================

  //  auto fd = io::open<ascii>("result.dat");
  balls = storage::handles<config::ball>(data, 0);
  auto ball1 = (balls | views::take(1)).front();
  auto ball2 = (balls | views::take(2)).back();

  // print("ball1:{}, ball2:{} ball1.q()={}, ball2.q()={}\n", ball1.get(),
  // ball2.get(), ball1.q(), ball2.q());
  //  fix this for constant fext
  simulation.one_step_integrator().compute_iteration_matrix(
      simulation.current_step());

  auto out = fmt::output_file("result-many.dat");

  out.print("{:.15e} {:.15e} {:.15e} {:.15e} {:.15e} {:.15e} {:.15e}\n",
            simulation.current_step() * simulation.time_step(),
            config::ball::q::at(ball1, simulation.current_step())(0),
            config::ball::q::at(ball2, simulation.current_step())(0),
            config::ball::velocity::at(ball1, simulation.current_step())(0),
            config::ball::velocity::at(ball2, simulation.current_step())(0),
            0., 0.);

  std::chrono::time_point<std::chrono::system_clock> start, end;
  start = std::chrono::system_clock::now();
  while (simulation.has_next_event()) {
    auto ninvds = simulation.compute_one_step();
    /*
        double p01, p02, lambda1, lambda2;
        if (ninvds > 1) {
          p01 =
       get_vector(simulation.one_step_integrator().p0_vector_assembled(),
                           0)(0);
          p02 =
       get_vector(simulation.one_step_integrator().p0_vector_assembled(),
                           0)(1);

          lambda1 = get_vector(
              simulation.one_step_integrator().lambda_vector_assembled(),
       0)(0);

          lambda2 = get_vector(
              simulation.one_step_integrator().lambda_vector_assembled(),
       1)(0);
        }
        else if (ninvds == 1) {
          p01 =
       get_vector(simulation.one_step_integrator().p0_vector_assembled(),
                           0)(0);
          p02 = 0;

          lambda1 = get_vector(
              simulation.one_step_integrator().lambda_vector_assembled(),
       0)(0);

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
                  config::ball::q::at(ball1, simulation.current_step())(0),
                  config::ball::q::at(ball2, simulation.current_step())(0),
                  config::ball::velocity::at(ball1,
       simulation.current_step())(0), config::ball::velocity::at(ball2,
       simulation.current_step())(0), p01, p02, lambda1, lambda2); */
  }

  end = std::chrono::system_clock::now();
  int elapsed =
      std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
          .count();
  print("RESULTS : -- {} --\n", elapsed);

  //  io::close(fd);
}
