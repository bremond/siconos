#include "siconos_data.hpp"
#include "siconos.hpp"

#include <tuple>
#include <fmt/core.h>
#include <fmt/ranges.h>
using fmt::print;


using namespace siconos;

int main()
{
  using formulation = lagrangian<param>;
  using ball = formulation::dynamical_system;
  using osnspb = one_step_nonsmooth_problem<lcp>;
  using osi = one_step_integrator<formulation>::moreau_jean;
  using relation = formulation::relation;
  using nslaw = nonsmooth_law::newton_impact;
  using interaction = interaction<nslaw, relation>;
  using td = time_discretization<param>;
  using simulation = time_stepping<td, osi, osnspb>;
  using siconos::get;
  auto data = siconos::make_data<env, simulation, interaction>();

  add_item<simulation>(0, data);
  get_memory<simulation::time_discretization::current_time_step>(data)[0][0] = 0;
  print("---\n");
  for_each(
    [](auto& a)
    {
      print("->{}\n", a.size());
    },
    data._collections);

  print("---\n");

  auto ds0 = add_vertex_item<ball>(0, data);
  auto& v0 = get<ball::velocity>(ds0, 0, data);

  v0 = { 1., 2., 3.};

  auto& velocityp = get<ball::velocity>(ds0, 0, data);
  print("--->{}\n", velocityp);

  for_each([](auto& a) { print("{:d}\n", a.size()); }, data._collections);
  print("{0}\n", data._collections);

  for_each([](auto& a) { print("{0}\n", a); }, data._collections);

  auto ds1 = add_vertex_item<ball>(0, data);

  get<ball::q>(ds1, 0, data) = { 1.,1., 1.};
  auto ds2 = add_vertex_item<ball>(0, data);
  get<ball::q>(ds2, 0, data) = { 9.,9., 9.};
  print("---\n");
  for_each([](auto& a) { print("{0}\n", a); }, data._collections);

  remove_vertex_item(ds1, 0, data);
  print("---\n");
  for_each([](auto& a) { print("{0}\n", a); }, data._collections);

  print("{}", get<ball::mass_matrix>(ds0, 0, data));

  auto m = get<ball::mass_matrix>(ds0, 0, data);

  print("---\n");

  print("{}", m);

  print("---\n");

  xget<ball::fext>(ds0, data) = { 10., 0., 0.};
  print("{}\n", data._collections);

  data(xget<ball::fext>)(ds0) = { 2.,1.,0.};
  print("{}\n", data(xget<ball::fext>)(ds0));

  add_item<nslaw>(0, data);

  auto& e = siconos::get_memory<nslaw::e>(data);
  e[0][0] = 0.7;

  print("{}\n", siconos::get_memory<nslaw::e>(data));
}

