#include "siconos_environment.hpp"
#include "siconos.hpp"

#include <tuple>
#include <fmt/core.h>
#include <fmt/ranges.h>
#include <functional>
using fmt::print;


using namespace siconos;
using env = siconos::standard_environment;

int main()
{
  using formulation = lagrangian<linear, time_invariant, degrees_of_freedom<3>>;
  using ball = formulation::dynamical_system;
  using osnspb = one_step_nonsmooth_problem<lcp>;
  using osi = one_step_integrator<formulation>::euler;
  using relation = formulation::relation;
  using nslaw = nonsmooth_law::newton_impact;
  using interaction = interaction<nslaw, relation>;
  using td = time_discretization<>;
  using topo = topology;
  using simulation = time_stepping<td, osi, osnspb, topo>;
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
  static_assert(match::attribute<nslaw::e>);

  static_assert(match::attribute_of<nslaw::e, nslaw>);
  static_assert(match::attribute_of<ball::velocity, ball>);
  static_assert(match::attribute_of<td::step, td>);

  static_assert(std::is_same_v<td, decltype(item_attribute<td::step>(all_items(simulation{})))>);

  static_assert(std::is_same_v<typename siconos::traits::config<env, typename siconos::time_discretization<>::step>::type,
                typename env::indice>);

  static_assert(filter<hold<decltype([]<typename T>(T){ return std::floating_point<T>; })>>(std::tuple<char,int,double>{}) == std::tuple<double>{});
//  static_assert(std::is_same_v<decltype(all_items_of_kind<graph_item>(simulation{})),
//                gather<ball, interaction>>);

//  static_assert(memory_size<typename td::step, typename osi::keeps> == 1);

//  static_assert(std::is_same_v<
//                decltype(all_keeps(simulation{})),
//                gather<keep<ball::q,2>, keep<ball::velocity,2>>>);

//  static_assert(memory_size<typename ball::q, typename osi::keeps> == 2);

  static_assert(std::array{1,2,3}-std::array{1,2,3}==std::array{0,0,0});
  static_assert(std::array{1,2,3}+std::array{1,2,3}==std::array{2,4,6});
  static_assert(std::array{std::array{1,2,3},std::array{4,5,6}}-std::array{std::array{1,2,3},std::array{4,5,6}}==std::array{std::array{0,0,0},std::array{0,0,0}});


  auto v1 = std::vector{1,2,3};
  auto v2 = std::vector{1,2,3};
  assert (v1+v2 == (std::vector{2,4,6}));

  auto v3 = std::vector{std::array{1,2,3},std::array{4,5,6}};
  auto v4 = std::vector{std::array{1,2,3},std::array{4,5,6}};
  assert(v3-v4 == (std::vector{std::array{0,0,0},std::array{0,0,0}}));
//  std::cout << boost::hana::experimental::type_name<decltype(attributes(interaction{}))>().c_str() << std::endl;
//  std::cout << boost::hana::experimental::type_name<gather<nslaw, relation>>().c_str() << std::endl;

//  std::cout << boost::hana::experimental::type_name<decltype(all_attributes(interaction{}))>().c_str() << std::endl;

  static_assert(std::is_same_v<decltype(attributes(interaction{})),
                gather<interaction::nonsmooth_law, interaction::relation>>);

  static_assert(std::is_same_v<decltype(all_attributes(interaction{})),
                gather<interaction::nonsmooth_law, interaction::relation,
                nslaw::e, relation::h_matrix>>);

  }

  auto data = siconos::make_storage<env, simulation,
                                    bounded_collection<relation, 3>,
                                    unbounded_collection<ball>,
                                    unbounded_collection<interaction>,
                                    with_kinds<
                                      keep<ball::q, 10>,
                                      diagonal<ball::mass_matrix>>>();
//   add<simulation>(data);
  print("v={}\n", siconos::get_memory<ball::velocity>(data));
  print("---\n");
  ground::for_each(data,
    [](auto& pair)
    {
      using value_t = std::decay_t<decltype(ground::second(pair))>;

      if constexpr (match::size<value_t>)
      {
        print("->{}\n", ground::second(pair).size());
      }
    });

   print("---\n");


   auto ds0 = fixed_add<ball>(data);

   ball::q::at(ds0) = { 0., 0., 1.};

   auto& v0 = ball::velocity::at(ds0);

   v0 = { 1., 2., 3.};
   print("v={}\n", siconos::get_memory<ball::velocity>(data));
   auto& velocityp = fix_map(get<ball::velocity>)(ds0);
   print("--->{}\n", velocityp);

   print("---\n");
   ground::for_each(data,
    [](auto& pair)
    {
      using value_t = std::decay_t<decltype(ground::second(pair))>;

      if constexpr (match::size<value_t>)
      {
        print("->{}\n", ground::second(pair).size());
      }
    });
   print("---\n");

   auto ds1 = fix(add<ball>)(data);
   ball::q::at(ds1) = { 1.,1., 1.};

   auto ds2 = fix(add<ball>)(data);
   ball::q::at(ds2) = { 9.,9., 9.};
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

   ball::fext::at(ds0) = { 10., 0., 0.};
//   print("{}\n", data._collections);

   ball::fext::at(ds0) = { 2.,1.,0.};
//   print("{}\n", data(get<ball::fext>)(ds0));

   auto inter1 = fix(add<interaction>)(data);
   auto nslaw1 = add<nslaw>(data);

   // siconos::set<interaction::nonsmooth_law>(inter1);
   interaction::nonsmooth_law::at(inter1) = nslaw1;
   fix_map(get<interaction::nonsmooth_law>)(inter1) = nslaw1;

   auto& e = siconos::get_memory<nslaw::e>(data);
   e[0][nslaw1.get()] = 0.7;

   print("{}\n", siconos::get_memory<nslaw::e>(data));

   auto inter2 = fix(add<interaction>)(data);
   auto nslaw2 = fix(add<nslaw>)(data);

   interaction::nonsmooth_law::at(inter2) = nslaw2;

   nslaw::e::internal_get(interaction::nonsmooth_law::at(inter2), data) = 0.3;

   auto ball1 = add<ball>(data);
   auto ball2 = add<ball>(data);
   interaction::nonsmooth_law::at(inter2);
   print("e={}\n", nslaw::e::at(nslaw2));

   using data_t = std::decay_t<decltype(data)>;

   print("memory_size={}\n", (memory_size<ball::q, decltype(all_kinds_of<some::keep>(data))>));
   print("q={}\n", siconos::get_memory<ball::q>(data));
   print("v={}\n", siconos::get_memory<ball::velocity>(data));

   auto ustore1 = unit_storage<env, some::scalar, some::vector<3>, some::graph>::type {};

//   ground::find(ustore1, boost::hana::type_c<some::scalar>) = 1.0;
//   ground::find(ustore1, boost::hana::type_c<some::vector<3>>).value() = { 1.0, 2.0, 3.0 };

   ground::get<some::scalar>(ustore1) = 2.0;
   ground::get<some::vector<3>>(ustore1) = std::array {1.0, 2.0, 3.0};

   print("ustore1 scalar={}\n", ground::get<some::scalar>(ustore1));

   print("ustore1 vector={}\n", ground::get<some::vector<3>>(ustore1));

   auto ustore2 = typename item_storage<env, simulation, ball>::type {};

   ground::get<ball::q>(ustore2) = { 2.0, 3.0, 4.0 };
   print("ustore2 ball::q ={}\n", ground::get<ball::q>(ustore2));

   //ground::transform(iget<ball>(ustore2), [&ustore2]<typename A>(A){ return iget<A>(ustore2); });

   auto ustore5 = make_storage<env, simulation,
                               unbounded_collection<ball>,
                               with_kinds<diagonal<ball::mass_matrix>>>();

   ground::get<ball::q>(ustore5)[0].push_back({7.,8.,9.});
   print("ustore5 ball::q ={}\n", ground::get<ball::q>(ustore5));

   auto some_iball=add<ball>(ustore5);

   get<ball::velocity>(some_iball, ustore5) = {3., 3., 3.};
   print("ustore5 ball::velocity = {}\n", get<ball::velocity>(some_iball, ustore5));

   // ok, compilation failure:
   // get<nslaw::e>(0, some_iball, ustore5);

   // TD: remove, nsds::link


   get<ball::fext>(some_iball, ustore5) = { 0., 0., 9.81 };

   auto the_td = add<td>(ustore5);
   get<td::step>(the_td, ustore5) = 0;
   get<td::h>(the_td, ustore5) = 0.01;

   auto m = 1.0;
   auto R = 1.0;

   print("mass matrix: {}", get<ball::mass_matrix>(some_iball, ustore5));
   get<ball::mass_matrix>(some_iball, ustore5) = { m, m , 2.5*m*R*R };
   print("mass matrix: {}", get<ball::mass_matrix>(some_iball, ustore5));

   simulation::compute_one_step(ustore5);
   print("current step : {}\n", simulation::current_step(ustore5));
   print("v={}\n", siconos::get_memory<ball::velocity>(ustore5));
   simulation::compute_one_step(ustore5);
   print("current step : {}\n", simulation::current_step(ustore5));
   print("v={}\n", siconos::get_memory<ball::velocity>(ustore5));
   simulation::compute_one_step(ustore5);
   print("current step : {}\n", simulation::current_step(ustore5));
   print("v={}\n", siconos::get_memory<ball::velocity>(ustore5));

   auto a = get<ball::velocity>(simulation::current_step(ustore5), some_iball, ustore5);
   print("ustore5 ball::velocity = {}\n", a);
}
