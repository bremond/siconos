#include "siconos_environment.hpp"
#include "siconos.hpp"

#include <tuple>
#include <fmt/core.h>
#include <fmt/ranges.h>

using fmt::print;


using namespace siconos;
using env = siconos::standard_environment;

int main()
{
  using formulation = lagrangian<linear, time_invariant, degrees_of_freedom<3>>;
  using ball = formulation::dynamical_system;
  using osnspb = one_step_nonsmooth_problem<lcp>;
  using osi = one_step_integrator<formulation>::moreau_jean;
  using relation = formulation::relation;
  using nslaw = nonsmooth_law::newton_impact;
  using interaction = interaction<nslaw, relation>;
  using td = time_discretization<>;
  using topo = topology;
  using simulation = time_stepping<td, osi, osnspb, topo, ball, interaction>;
  using siconos::get;

  {
  static_assert(must::contains<int,std::tuple<int,float,char>>);
  static_assert(!must::contains<double,std::tuple<int,float,char>>);
  static_assert(std::is_same_v<decltype(
                  cons(int{}, std::tuple<char, float>{})),
                std::tuple<int,char,float>>);

  static_assert(std::is_same_v<
                decltype(append(
                           gather<int>{}, gather<float>{})),
                gather<int, float>>);

  static_assert(transform([]<typename T>(T) { return int{}; }, std::tuple<char,float,double>{}) == std::tuple<int,int,int>{});
  static_assert(std::is_same_v<decltype(all_items(nslaw{})), gather<siconos::nonsmooth_law::newton_impact>>);
  static_assert(
    std::is_same_v<decltype([]()
    {
      return transform([]<typename T>(T) { return typename T::type{}; },
                       filter<hold<decltype([]<typename T>(T) { return match::item_ref<T>; })>>(typename interaction::attributes{}));
    }()), gather<nslaw, relation>>);

  static_assert(std::is_same_v<decltype(all_items(interaction{})),
                std::tuple<siconos::interaction<siconos::nonsmooth_law::newton_impact, siconos::lagrangian<linear, time_invariant, degrees_of_freedom<3>>::relation>, siconos::nonsmooth_law::newton_impact, siconos::lagrangian<linear, time_invariant, degrees_of_freedom<3>>::relation>>);

  static_assert(must::contains<osnspb, decltype(all_items(simulation{}))>);

  static_assert(match::item<ball>);
  static_assert(match::vertex_item<ball>);
  static_assert(match::attribute<nslaw::e>);

  static_assert(match::attribute_of<nslaw::e, nslaw>);
  static_assert(match::attribute_of<ball::velocity, ball>);
  static_assert(match::attribute_of<td::step, td>);

  static_assert(std::is_same_v<td, decltype(item_attribute<td::step>(all_items(simulation{})))>);

  static_assert(std::is_same_v<typename siconos::traits::config<env, typename siconos::time_discretization<>::step>::type,
                typename env::indice>);
  static_assert(attribute_storage_t<
                env,
                td::step,
                decltype(all_items(simulation{}))>{} == typename place_holder<typename env::indice>::type{});
  static_assert(std::is_same_v<attribute_storage_t<
                env,
                td::step,
                decltype(all_items(simulation{}))>,
                typename place_holder<typename env::indice>::type>);

  static_assert(filter<hold<decltype([]<typename T>(T){ return std::floating_point<T>; })>>(std::tuple<char,int,double>{}) == std::tuple<double>{});
  static_assert(std::is_same_v<decltype(all_graph_items(simulation{})),
                gather<ball, interaction>>);

  static_assert(memory_size<typename td::step, typename osi::keeps> == 1);

  static_assert(std::is_same_v<
                decltype(all_keeps(simulation{})),
                gather<keep<ball::q,2>, keep<ball::velocity,2>>>);

  static_assert(memory_size<typename ball::q, typename osi::keeps> == 2);
  }
  auto data = siconos::make_storage<env, simulation>();
//   add<simulation>(data);

  print("---\n");
  ground::for_each(data._collections,
    [](auto& a)
    {
      print("->{}\n", a.size());
    });

   print("---\n");


   auto ds0 = fixed_add<ball>(data);

   ball::q::get(ds0) = { 0., 0., 1.};

   auto& v0 = ball::velocity::get(ds0);

   v0 = { 1., 2., 3.};

   auto& velocityp = fix_map(get<ball::velocity>)(ds0);
   print("--->{}\n", velocityp);

   ground::for_each(data._collections, [](auto& a) { print("{:d}\n", a.size()); });
//   print("{0}\n", data._collections);

//   for_each([](auto& a) { print("{0}\n", a); }, data._collections);

   auto ds1 = fix(add<ball>)(data);
   ball::q::get(ds1) = { 1.,1., 1.};

   auto ds2 = fix(add<ball>)(data);
   ball::q::get(ds2) = { 9.,9., 9.};
//   print("---\n");
//   for_each([](auto& a) { print("{0}\n", a); }, data._collections);

//   remove_vertex_item(ds1, 0);
//   print("---\n");
//   for_each([](auto& a) { print("{0}\n", a); }, data._collections);

   print("-----{}------", fix_map(get<ball::mass_matrix>)(ds1));
   print("{}", fix_map(get<ball::mass_matrix>)(ds0));

//   auto m = get<ball::mass_matrix>(ds0, data);

//   print("---\n");

//   print("{}", m);

//   print("---\n");

   ball::fext::get(ds0) = { 10., 0., 0.};
//   print("{}\n", data._collections);

   ball::fext::get(ds0) = { 2.,1.,0.};
//   print("{}\n", data(get<ball::fext>)(ds0));

   auto inter1 = fix(add<interaction>)(data);
   auto nslaw1 = add<nslaw>(data);

   // siconos::set<interaction::nonsmooth_law>(inter1);
   interaction::nonsmooth_law::get(inter1) = nslaw1;
   fix_map(get<interaction::nonsmooth_law>)(inter1) = nslaw1;

   auto& e = siconos::get_memory<nslaw::e>(data);
   e[0][nslaw1.get()] = 0.7;

   print("{}\n", siconos::get_memory<nslaw::e>(data));

   auto inter2 = fix(add<interaction>)(data);
   auto nslaw2 = fix(add<nslaw>)(data);

   interaction::nonsmooth_law::get(inter2) = nslaw2;

   nslaw::e::internal_get(interaction::nonsmooth_law::get(inter2), data) = 0.3;

   auto ball1 = add<ball>(data);
   auto ball2 = add<ball>(data);
   interaction::nonsmooth_law::get(inter2);
   print("e={}\n", nslaw::e::get(nslaw2));

   using data_t = std::decay_t<decltype(data)>;

   print("memory_size={}\n", (memory_size<ball::q, typename data_t::machine::osi::keeps> ));
   print("q={}\n", siconos::get_memory<ball::q>(data));
   print("v={}\n", siconos::get_memory<ball::velocity>(data));

   auto ustore1 = unit_storage<env, some::scalar, some::vector<3>, some::graph>::type {};

//   ground::find(ustore1, boost::hana::type_c<some::scalar>) = 1.0;
//   ground::find(ustore1, boost::hana::type_c<some::vector<3>>).value() = { 1.0, 2.0, 3.0 };

   ground::get<some::scalar>(ustore1) = 2.0;
   ground::get<some::vector<3>>(ustore1) = std::array {1.0, 2.0, 3.0};

   print("ustore1 scalar={}\n", ground::get<some::scalar>(ustore1));

   print("ustore1 vector={}\n", ground::get<some::vector<3>>(ustore1));


   auto ustore2 = std::apply([]<typename ...Args>(Args...)
                             { return typename unit_storage<env, Args...> ::type {};}, all_attributes(simulation{}));

   ground::get<ball::q>(ustore2) = { 2.0, 3.0, 4.0 };

   print("ustore2 ball::q ={}\n", ground::get<ball::q>(ustore2));

   auto ustore3 = ground::map_vtransform([](auto&&s) { return std::vector {s}; })(ustore2);

   ground::get<ball::q>(ustore3).push_back({4.,5.,6.});

   print("ustore3 ball::q ={}\n", ground::get<ball::q>(ustore3));

}
