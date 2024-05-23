#include <siconos/numerics/Friction_cst.h>

#include "siconos/collision/space_filter.hpp"
#include "siconos/siconos.hpp"
#include "siconos/utils/print.hpp"

#include <bits/stdc++.h>

namespace siconos::config {

using disk = model::lagrangian_ds;
using fc2d = simul::nonsmooth_problem<FrictionContactProblem>;
using osnspb = simul::one_step_nonsmooth_problem<fc2d>;
using nslaw = model::newton_impact_friction;
using disk_shape = collision::shape::disk;
using diskdisk_r = collision::diskdisk_r;
using disksegment_r = collision::disksegment_r;
using interaction = simul::interaction<nslaw, diskdisk_r, disksegment_r>;
using osi = simul::one_step_integrator<disk, interaction>::moreau_jean;
using td = simul::time_discretization<>;
using topo = simul::topology<disk, interaction>;
using simulation = simul::time_stepping<td, osi, osnspb, topo>;
using pointd = collision::point<disk>;
using pointl = collision::point<collision::shape::segment>;
using neighborhood = collision::neighborhood<pointd, pointl>;
using space_filter = collision::space_filter<topo, neighborhood>;

using params = map<iparam<"dof", 3>>;
}  // namespace siconos::config

int main(int argc, char* argv[])
{
  using namespace siconos;
  using storage::pattern::wrap;
  using namespace storage;

  auto data = storage::make<
      standard_environment<config::params>, config::simulation,
      wrap<some::unbounded_collection, config::disk>, config::disk_shape,
      config::diskdisk_r,
      wrap<some::unbounded_collection, config::disksegment_r>,
      wrap<some::unbounded_collection, config::pointl>,
      wrap<some::unbounded_collection, config::pointd>,
      wrap<some::unbounded_collection, config::interaction>,
      config::space_filter,
      storage::with_properties<
          storage::attached<config::disk, storage::pattern::symbol<"shape">,
                            storage::some::item_ref<config::disk_shape>>,
          storage::time_invariant<storage::attr_t<config::disk, "fext">>,
          storage::diagonal<storage::attr_t<config::disk, "mass_matrix">>,
          storage::unbounded_diagonal<
              storage::attr_t<config::osi, "mass_matrix_assembled">>>>();

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
//  auto d2 = storage::add<config::disk>(data);

  d1.q() = {0, position_init, 0};
  d1.velocity() = {0, velocity_init, 0};
  d1.mass_matrix().diagonal() << m, m, 2. / 5. * m * radius * radius;

  // d2.q() = {2 * position_init, 0.1, 0};
  // d2.velocity() = {velocity_init, 0, 0};
  // d2.mass_matrix().diagonal() << m, m, 2. / 5. * m * radius * radius;

  storage::handle(data, prop<"shape">(d1)).radius() = radius;
//  storage::handle(data, prop<"shape">(d2)).radius() = 2;

  // -- Set external forces (weight) --
  d1.fext() = {0., -m * g, 0.};
//  d2.fext() = {-m * g, 0., 0.};

  // ------------------
  // -- The relation --
  // ------------------

  // -- nslaw --
  double e = 0.9;
  auto nslaw = storage::add<config::nslaw>(data);
  nslaw.e() = e;

  auto fc2d = storage::add<config::fc2d>(data);
  fc2d.create();
  fc2d.instance()->dimension = 2;
//  fc2d.instance()->mu = 0.1;

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
  so.create(SICONOS_FRICTION_2D_LEMKE);
  osnspb.options() = so;

  auto ngbh = storage::add<config::neighborhood>(data);

  ngbh.create(2.);  // radius

  auto diskdisk_r = storage::add<config::diskdisk_r>(data);
  auto ground_r = storage::add<config::disksegment_r>(data);
  auto segment = storage::handle(data, ground_r.segment());
  segment.x1() = 0.;
  segment.y1() = 0.;
  segment.x2() = 10.;
  segment.y2() = 0.;
  segment.initialize();
  segment.maxpoints() = 100;

  auto spacef = storage::add<config::space_filter>(data);
  spacef.neighborhood() = ngbh;
  spacef.diskdisk_r() = diskdisk_r;
  spacef.nslaw() = nslaw;
  spacef.insert_disksegment_r(ground_r);
  spacef.make_points();
  ngbh.add_point_sets(0);
  // =========================== End of model definition
  // ===========================
  // ================================= Computation
  // =================================

  //  auto fd = io::open<ascii>("result.dat");

  // fix this for constant fext
  simulation.initialize();

//  auto out = fmt::output_file("result.dat");
  std::ofstream cout("result.dat");

  // https://stackoverflow.com/questions/72767354/how-to-flush-fmt-output-in-debug-mode
  cout << fmt::format("{:.15e} {:.15e} {:.15e} {:.15e} {:.15e}\n",
            simulation.current_step() * simulation.time_step(),
            storage::attr<"q">(d1, simulation.current_step())(1),
            storage::attr<"velocity">(d1, simulation.current_step())(1), 0.,
                          0.) << std::flush;

  while (simulation.has_next_event()) {
    ngbh.update(0);
    ngbh.search();
    spacef.update_index_set0();

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

    cout << fmt::format("{:.15e} {:.15e} {:.15e} {:.15e} {:.15e}\n",
              simulation.current_step() * simulation.time_step(),
              storage::attr<"q">(d1, simulation.current_step())(1),
              storage::attr<"velocity">(d1, simulation.current_step())(1), p0,
                            lambda) << std::flush;
  }
  //  io::close(fd);
}
