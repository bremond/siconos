#include <siconos/siconos.hpp>

using namespace siconos;

int main(int argc, char* argv[])
{
  using formulation = lagrangian<linear, time_invariant, dof<3>>;
  using osi = one_step_integrator<formulation>::moreau_jean;
  using osnspb = one_step_nonsmooth_problem<lcp>;
  using ball = formulation::dynamical_system;
  using nslaw = nonsmooth_law::newton_impact;
  using relation = relation<formulation>;
  using interaction = interaction<nslaw, relation>;
  using simulation = time_stepping<osi, osnspb>;

  auto data = make_data<std_env, simulation, interaction>;

  //unsigned int nDof = 3;           // degrees of freedom for the ball
  double t0 = 0;                   // initial computation time
  double T = 10;                  // final computation time
  double h = 0.005;                // time step
  double position_init = 1.0;      // initial position for lowest bead.
  double velocity_init = 0.0;      // initial velocity for lowest bead.
  double theta = 0.5;              // theta for MoreauJeanOSI integrator
  double R = 0.1; // Ball radius
  double m = 1.; // Ball mass
  double g = 9.81; // Gravity

  print("====> Model loading ...\n");

  // -- The dynamical_system --
  auto the_ball = add<ball>(data);

  get<ball::q>(the_ball, data) = {position_init, 0, 0};
  // set<ball::q>(the_ball, data, {position_init, 0, 0});
  // data(get<ball::q>) = {position_init, 0, 0};
  get<ball::velocity>(the_ball, data) = {position_init, 0, 0};
  get<ball::mass_matrix>(the_ball, data) =
    {m, 0, 0
     0, m, 0,
     0, 0, 2./5.*m*R*R};;

  // -- Set external forces (weight) --
  get<ball::fext>(the_ball, data) = {-m*g, 0., 0.};

  // --------------------
  // --- Interactions ---
  // --------------------

  // Lagrangian relation
  auto the_relation = add<relation>(data);
  get<relation::h_matrix>(the_relation, data) = {1.0, 0., 0.};

  // Interaction ball-floor
  auto the_interaction = add<interaction>(data);

  // -- nslaw --
  get<nslaw::e>(the_interaction, data) = e;

  get<interaction::relation>(the_interaction, data) = { the_relation };
  get<interaction::link>(the_interaction, data) = { the_ball };

  // ------------------
  // --- Simulation ---
  // ------------------
  get<osi::theta>(data) = theta;

  get<time_discretization::t0>(data) = t0;
  get<time_discretization::h>(data) = h;

  // =========================== End of model definition ===========================
  // ================================= Computation =================================

  auto fd = io::open<ascii>("result.dat");

  while (simulation::has_next_event(data))
  {
    simulation::compute_one_step(data);

    io::write(get<ball::q, 0>,
              get<ball::v, 0>,
              get<ball::p, 1>
              get<interaction::lambda, 1>,
              { the_ball }, fd, data);

  }

  io::close(fd);
}
