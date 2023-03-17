#include "siconos_environment.hpp"
#include "siconos.hpp"

#include <fmt/core.h>
#include <fmt/ranges.h>

using namespace siconos;

using fmt::print;
int main(int argc, char* argv[])
{
  using formulation = lagrangian<linear, time_invariant, degrees_of_freedom<3>>;
  using osnspb = one_step_nonsmooth_problem<lcp>;
  using ball = formulation::dynamical_system;
  using nslaw = nonsmooth_law::newton_impact;
  using relation = formulation::relation;
  using interaction = interaction<nslaw, formulation, 1>;
  using osi = one_step_integrator<formulation, interaction>::moreau_jean;
  using td = time_discretization<>;
  using topo = topology<formulation, interaction>;  // topology<ball, interaction>
  using simulation = time_stepping<td, osi, osnspb, topo>;
  using siconos::get;

  auto data = make_storage<standard_environment, simulation, ball, interaction>();

  //unsigned int nDof = 3;         // degrees of freedom for the ball
  double t0 = 0;                   // initial computation time
  double T = 10;                   // final computation time
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
  get<ball::velocity>(the_ball, data) = {velocity_init, 0, 0};
  get<ball::mass_matrix>(the_ball, data) =
    {m, 0, 0,
     0, m, 0,
     0, 0, 2./5.*m*R*R};

  // -- Set external forces (weight) --
  get<ball::fext>(the_ball, data) = {-m*g, 0., 0.};

  // --------------------
  // --- Interactions ---
  // --------------------

  // -- Lagrangian relation --
  auto the_relation = add<relation>(data);
  get<relation::h_matrix>(the_relation, data) = {1.0, 0., 0.};

  // -- nslaw --
  double e = 0.9;
  auto the_nslaw = add<nslaw>(data);
  get<nslaw::e>(the_nslaw, data) = e;

  // Interaction ball-floor
  auto the_interaction = add<interaction>(data);

  get<interaction::relation>(the_interaction, data) =  the_relation ;
  get<interaction::nonsmooth_law>(the_interaction, data) = the_nslaw ;
  //get<interaction::link>(the_interaction, data) =  the_ball ;

  // ------------------
  // --- Simulation ---
  // ------------------
  auto the_osi = add<osi>(data);
  auto the_td = add<td>(data);
  auto the_osnspb = add<osnspb>(data);
  auto the_simulation = add<simulation>(data);

  // get<simulation::osi>(the_simulation, data) = the_osi;
  // get<simulation::td>(the_simulation, data) = the_td;
  // get<simulation::osnspb>(the_simulation, data) = the_osnspb;

  get<osi::theta>(the_osi, data) = theta;
  get<td::t0>(the_td, data) = t0;
  get<td::h>(the_td, data) = h;

  // =========================== End of model definition ===========================
  // ================================= Computation =================================

#if defined(NONONO)
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

#endif
}

