#include "siconos/utils/check.hpp"
#include "siconos/utils/environment.hpp"
#include "siconos/siconos.hpp"
#include "siconos/utils/ground.hpp"
#include "siconos/utils/pattern.hpp"
#include "siconos/algebra/numerics.hpp"
#include "siconos/model/lagrangian_ds.hpp"
#include "siconos/model/lagrangian_r.hpp"
#include "siconos/model/nslaws.hpp"
#include "siconos/simul/topology.hpp"
#include "siconos/simul/interaction.hpp"
#include "siconos/simul/one_step_nonsmooth_problem.hpp"
#include "siconos/simul/one_step_integrator.hpp"
#include "siconos/simul/time_discretization.hpp"
#include "siconos/simul/time_stepping.hpp"
#include <charconv>
#include <tuple>
#include <fmt/core.h>
#include <fmt/ranges.h>
#include <functional>
#include <typeinfo>
#if defined(__clang__)
#include <boost/hana/experimental/type_name.hpp>
#endif
using fmt::print;

using namespace siconos;
using env = siconos::standard_environment;


using namespace boost::hana::literals;
int main()
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

  struct zz {};
//  {
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

  //ground::type_trace<decltype(all_items(interaction{}))>();
  static_assert(std::is_same_v<decltype(all_items(interaction{})),
                std::tuple<interaction, nslaw, relation>>);

  //static_assert(must::contains<osnspb, decltype(all_items(simulation{}))>);

  static_assert(match::item<ball>);
  static_assert(match::attribute<nslaw::e>);

  static_assert(match::attribute_of<nslaw::e, nslaw>);
  static_assert(match::attribute_of<ball::velocity, ball>);
  static_assert(match::attribute_of<td::step, td>);

  //static_assert(std::is_same_v<td, decltype(item_attribute<td::step>(all_items(simulation{})))>);

  static_assert(std::is_same_v<typename siconos::traits::config<env>::convert<typename siconos::time_discretization<>::step>::type,
                typename env::indice>);

  static_assert(std::is_same_v<typename siconos::traits::config<env>::convert<some::scalar>::type,
                typename env::scalar>);

  static_assert(match::attached_storage<attached_storage<ball, zz, some::scalar>, ball>);
  static_assert(match::attached_storage<attached_storage<ball, zz, some::scalar>,
                wrap<some::unbounded_collection, ball>>);
  static_assert(match::attached_storage<attached_storage<ball, symbol<"Z">, some::scalar>, wrap<some::unbounded_collection, ball>>);
  static_assert(!match::attached_storage<attached_storage<ball, zz, some::scalar>, nslaw>);


  static_assert(match::unbounded_storage<some::unbounded_collection<some::scalar>>);
  static_assert(match::bounded_storage<some::bounded_collection<some::scalar, some::indice_value<1>>>);
  static_assert(!match::unbounded_storage<some::bounded_collection<some::scalar, some::indice_value<1>>>);
  static_assert(!match::bounded_storage<some::unbounded_collection<some::scalar>>);
  static_assert(match::unbounded_storage<some::unbounded_diagonal_matrix<some::scalar>>);

  static_assert(match::wrap<wrap<some::unbounded_diagonal_matrix, ball>>);

  static_assert(traits::translatable<int, env>);

  std::cout << ground::type_name<traits::config<env>::convert<int>::type>();
  static_assert(std::is_same_v<traits::config<env>::convert<int>::type, int>);

  static_assert(traits::translatable<siconos::some::unbounded_collection<char[3]>, env>);
  static_assert(filter<hold<decltype([]<typename T>(T){ return std::floating_point<T>; })>>(std::tuple<char,int,double>{}) == std::tuple<double>{});

  static_assert(ground::filter(std::tuple<char, int, double>{},
                               ground::compose(
                                   ground::trait<std::is_floating_point>,
                                   ground::typeid_)) == std::tuple<double>{});

  static_assert(
      std::is_same_v <
          decltype(ground::filter(
              std::tuple<nslaw::e, ball::velocity>{},
              ground::derive_from<some::scalar>)),
      std::tuple<nslaw::e>>);
  //  static_assert(std::is_same_v<decltype(all_items_of_kind<graph_item>(simulation{})),
  //                gather<ball, interaction>>);

  //  static_assert(memory_size<typename td::step, typename osi::keeps> == 1);

  //  static_assert(std::is_same_v<
  //                decltype(all_keeps(simulation{})),
  //                gather<keep<ball::q,2>, keep<ball::velocity,2>>>);

  //  static_assert(memory_size<typename ball::q, typename osi::keeps> == 2);

//  static_assert(std::array{1,2,3}-std::array{1,2,3}==std::array{0,0,0});
//  static_assert(std::array{1,2,3}+std::array{1,2,3}==std::array{2,4,6});
//  static_assert(std::array{std::array{1,2,3},std::array{4,5,6}}-std::array{std::array{1,2,3},std::array{4,5,6}}==std::array{std::array{0,0,0},std::array{0,0,0}});


//  auto v1 = std::vector{1,2,3};
//  auto v2 = std::vector{1,2,3};
//  assert (v1+v2 == (std::vector{2,4,6}));

//  auto v3 = std::vector{std::array{1,2,3},std::array{4,5,6}};
//  auto v4 = std::vector{std::array{1,2,3},std::array{4,5,6}};
//  assert(v3-v4 == (std::vector{std::array{0,0,0},std::array{0,0,0}}));
//  std::cout << boost::hana::experimental::type_name<decltype(attributes(interaction{}))>().c_str() << std::endl;
//  std::cout << boost::hana::experimental::type_name<gather<nslaw, relation>>().c_str() << std::endl;

//  std::cout << boost::hana::experimental::type_name<decltype(all_attributes(interaction{}))>().c_str() << std::endl;

  static_assert(std::is_same_v<decltype(attributes(interaction{})),
                gather<some::indice_parameter<"dof">, some::indice_value<1>,interaction::nonsmooth_law, interaction::relation, interaction::h_matrix, interaction::lambda, interaction::y>>);

  static_assert(std::is_same_v<decltype(all_attributes(interaction{})),
                gather<some::indice_parameter<"dof">, some::indice_value<1>,interaction::nonsmooth_law, interaction::relation, interaction::h_matrix,
                interaction::lambda, interaction::y, nslaw::e>>);

  //}

   auto ddd = make_storage<env, bbb>();

   auto bob1 = add<bbb>(ddd);

   bbb::attr::at(bob1).reset(new aaa);

  auto data0 = siconos::make_storage<env, nslaw>();
  auto nslaw0 = add<nslaw>(data0);

  get<nslaw::e>(nslaw0, data0) = 0.6;
  get<nslaw::e>(nslaw0, data0) = 0.7;
  assert(get<nslaw::e>(nslaw0, data0) == 0.7);

  auto x = get<nslaw::e>(nslaw0, data0);
  x = 0.4;

  nslaw::e::at(nslaw0) = 0.8;

  assert (get<nslaw::e>(nslaw0, data0) == 0.8);

  auto data1 = siconos::make_storage<env, wrap<some::unbounded_collection, nslaw>>();

  auto nslaw1 = add<nslaw>(data1);

  nslaw::e::at(nslaw1) = 0.9;

  auto data2 = siconos::make_storage<env,
                                     wrap<some::unbounded_collection, nslaw>,
                                     with_properties<attached_storage<nslaw, zz, some::scalar>>>();

  auto nslaw2b = add<nslaw>(data2);

  nslaw2b.property(zz{}) = 2;

  auto data = siconos::make_storage<env, simulation,
                                    wrap<some::unbounded_collection, relation>,
                                    wrap<some::unbounded_collection, ball>,
                                    wrap<some::unbounded_collection, interaction>,
                                    with_properties<
                                      keep<ball::q, 10>,
                                      diagonal<ball::mass_matrix>,
                                      unbounded_diagonal<osi::mass_matrix_assembled>>>();
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


   auto ds0 = add<ball>(data);

//   ground::type_trace<decltype(ball::fext::at(ds0))>();
   ball::fext::at(ds0) = { 10., 0., 0.};
   ball::q::at(ds0) = { 0., 0., 1.};

   auto& v0 = ball::velocity::at(ds0);

   v0 = { 1., 2., 3.};
   print("v={}\n", siconos::get_memory<ball::velocity>(data));
   auto& velocityp = get<ball::velocity>(ds0, data);
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

   auto ds1 = add<ball>(data);
   ball::fext::at(ds0) = { 10., 0., 0.};
   ball::q::at(ds1) = { 1.,1., 1.};
   ball::fext::at(ds0) = { 10., 0., 0.};

   auto ds2 = add<ball>(data);
   ball::fext::at(ds0) = { 10., 0., 0.};
   ball::q::at(ds2) = { 9.,9., 9.};
   ball::fext::at(ds0) = { 10., 0., 0.};
//   print("---\n");
//   for_each([](auto& a) { print("{0}\n", a); }, data._collections);

//   remove_vertex_item(ds1, 0);
//   print("---\n");
//   for_each([](auto& a) { print("{0}\n", a); }, data._collections);

   //rint("-----{}------", get<ball::mass_matrix>(ds1, data));
   //print("{}", get<ball::mass_matrix>(ds0, data));

//   auto m = get<ball::mass_matrix>(ds0, data);

//   print("---\n");

//   print("{}", m);

//   print("---\n");

   ball::fext::at(ds0) = { 10., 0., 0.};
//   print("{}\n", data._collections);

   ball::fext::at(ds0) = { 2.,1.,0.};
//   print("{}\n", data(get<ball::fext>)(ds0));

   auto inter1 = add<interaction>(data);
   auto nslaw2 = add<nslaw>(data);

   auto dsimul = add<simulation>(data);

   // siconos::set<interaction::nonsmooth_law>(inter1);
   interaction::nonsmooth_law::at(inter1) = nslaw2;
   get<interaction::nonsmooth_law>(inter1, data) = nslaw2;

   auto& e = siconos::get_memory<nslaw::e>(data);
   e[0][nslaw2.get()] = 0.7;

   nslaw2.e() = 0.6;

   print("{}\n", siconos::get_memory<nslaw::e>(data));

   auto inter2 = add<interaction>(data);
   auto nslaw3 = add<nslaw>(data);

   interaction::nonsmooth_law::at(inter2) = nslaw3;
   nslaw::e::at(interaction::nonsmooth_law::at(inter2), data) = 0.3;

   auto htopo = dsimul.topology();
   auto hosi = dsimul.one_step_integrator();
   auto ball1 = add<ball>(data);
   auto ball2 = add<ball>(data);
   auto ball3 = add<ball>(data);
   print("e={}\n", nslaw::e::at(nslaw3));

   auto intera = htopo.link(ball1);
   auto interb = htopo.link(ball1, ball2);
   auto interc = htopo.link(ball1, ball3);
   auto interd = htopo.link(ball3);

   print("intera.nds={}\n", intera.property(symbol<"nds">{}));
   print("interb.nds={}\n", interb.property(symbol<"nds">{}));
   print("interc.nds={}\n", interc.property(symbol<"nds">{}));
   print("interd.nds={}\n", interd.property(symbol<"nds">{}));

   handle(htopo.dynamical_system_graphs()[0].bundle(ball1.property(symbol<"vd">{})), data).property(symbol<"index">{}) +=10;;

   ball2.property(symbol<"index">{}) +=1;
   print("ball1.index = {}\n", ball1.property(symbol<"index">{}));
   print("ball2.index = {}\n", ball2.property(symbol<"index">{}));

   print("htopo.index : {}\n", ground::get<attached_storage<ball, symbol<"index">, some::indice>>(data));
   print("make_index\n");
   htopo.make_index();
   print("htopo.index : {}\n", ground::get<attached_storage<ball, symbol<"index">, some::indice>>(data));

   hosi.assemble_h_matrix_for_involved_ds(0);
   hosi.assemble_mass_matrix_for_involved_ds(0);

   hosi.compute_w_matrix(0);
   hosi.compute_q_vector_assembled(0);

   auto so = add<numerics::solver_options>(data);
   so.create();

   auto xso = so.instance();
   print("solver options iSize = {}\n", xso->iSize);
   print("SICONOS_IPARAM_MAX_ITER = {}\n", SICONOS_IPARAM_MAX_ITER);
   print("iparam[SICONOS_IPARAM_MAX_ITER] = {}\n", xso->iparam[SICONOS_IPARAM_MAX_ITER]);

   so.instance()->iparam[SICONOS_IPARAM_MAX_ITER] = 100;

   auto lcp_1 = add<lcp>(data);
   lcp_1.create();

   auto osnspb_1 = dsimul.one_step_nonsmooth_problem();
   osnspb_1.problem() = lcp_1;
   osnspb_1.options() = so;

   print("memory_size={}\n", (memory_size<ball::q, decltype(all_properties_as<property::keep>(data))>));

   hosi.compute_q_vector_assembled(0);
   dsimul.solve_nonsmooth_problem<LinearComplementarityProblem>();
//   auto& q = siconos::get_memory<ball::q>(data);
//   print("q={}\n", q);
//   print("v={}\n", siconos::get_memory<ball::velocity>(data));

   auto ustore1 = unit_storage<env, some::scalar, some::vector<some::scalar, some::indice_value<3>>, some::graph<some::indice, some::indice>>::type {};

//   ground::find(ustore1, boost::hana::type_c<some::scalar>) = 1.0;
//   ground::find(ustore1, boost::hana::type_c<some::vector<3>>).value() = { 1.0, 2.0, 3.0 };

   ground::get<some::scalar>(ustore1) = 2.0;
   ground::get<some::vector<some::scalar, some::indice_value<3>>>(ustore1) =  {1.0, 2.0, 3.0};

   print("ustore1 scalar={}\n", ground::get<some::scalar>(ustore1));

   print("ustore1 vector={}\n", ground::get<some::vector<some::scalar, some::indice_value<3>>>(ustore1));

   auto ustore2 = typename item_storage<env, simulation, ball>::type {};

   ground::get<ball::q>(ustore2) = { 2.0, 3.0, 4.0 };
   print("ustore2 ball::q ={}\n", ground::get<ball::q>(ustore2));

   //ground::transform(iget<ball>(ustore2), [&ustore2]<typename A>(A){ return iget<A>(ustore2); });

   auto ustore5 = make_storage<env, simulation,
                               wrap<some::unbounded_collection, ball>,
                               wrap<some::unbounded_collection, interaction>,
                               with_properties<diagonal<ball::mass_matrix>,
                                               attached_storage<ball, zz , some::scalar>>>();


   //ground::get<ball::q>(ustore5)[0].push_back({7.,8.,9.});
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
//   print("mass matrix: {}", get<ball::mass_matrix>(some_iball, ustore5));
   auto ma = some_iball.mass_matrix();//get<ball::mass_matrix>(some_iball, ustore5);
//   ma.resize(3,3);
   ma.diagonal()[0] = m;
   ma.diagonal()[1] = m;
   ma.diagonal()[2] = 2.5*m*R*R;
//   print("mass matrix: {}", get<ball::mass_matrix>(some_iball, ustore5));

   auto simul = make_full_handle<simulation>(0, ustore5);

   simul.update_indexsets(0);
   simul.compute_one_step();

   print("current step : {}\n", simul.current_step());
   print("v={}\n", siconos::get_memory<ball::velocity>(ustore5));
   simul.compute_one_step();
   print("current step : {}\n", simul.current_step());
   print("v={}\n", siconos::get_memory<ball::velocity>(ustore5));
   simul.compute_one_step();
   print("current step : {}\n", simul.current_step());
   print("v={}\n", siconos::get_memory<ball::velocity>(ustore5));

   auto a = get<ball::velocity>(simul.current_step(), some_iball, ustore5);
   print("ustore5 ball::velocity = {}\n", a);

   auto& dsg0 = get<topo::dynamical_system_graphs>(htopo, ustore5)[0];

   auto dsgv = dsg0.add_vertex(some_iball);

   auto& ig0 = get<topo::interaction_graphs>(htopo, ustore5)[0];

   auto inter = add<interaction>(ustore5);

   auto [new_edge, ig_new_ve] = dsg0.add_edge(dsgv, dsgv, inter, ig0);
   auto rel = add<relation>(ustore5);

   auto yyv = ground::get<interaction::h_matrix>(ustore5);
   interaction::h_matrix::at(inter)(0, 0) = 1;

   interaction::relation::at(inter) = rel;
   auto& v = interaction::y::at(inter)[0];
   v = { 1 };

//   rel.compute_output(0., some_iball, inter, 1_c);

//   rel.compute_input(0., some_iball, inter, 1_c);
   auto hsim = make_full_handle<simulation>(0, data);
//   hsim.compute_output(1_c);

   ground::get<attached_storage<ball, zz, some::scalar>>(ustore5)[0][0] = 1.0;
   get<zz>(0, some_iball, ustore5) = 1.0;
   some_iball.property(zz{}) = 1.0;
//   for_each(ustore5,
   //           []<typename Key, typename Value>(Key&& key, Value&& value)
   //        {
//              print("{}\n", boost::hana::experimental::type_name<Key>().c_str());
   //       });


#ifdef __clang__
//   print("-->{}\n", boost::hana::experimental::type_name<ground::pair<ball, zz>>().c_str());
#endif
   auto dd = make_storage<env, wrap<some::unbounded_collection, ball>,
                          with_properties<diagonal<ball::mass_matrix>, keep<ball::q, 10>,
                                          attached_storage<ball, decltype("z"_s), some::scalar>,
                                          attached_storage<ball, symbol<"x">, some::scalar>>>();

  auto xball1 = add<ball>(dd);
  auto xball2 = add<ball>(dd);
  auto xball3 = add<ball>(dd);

  print("{},{}\n", xball1.get(), xball2.get());
  xball1.velocity() = {1, 1, 1};
  xball2.velocity() = {2, 2, 2};
  xball3.velocity() = {3, 3, 3};
  print("{}\n", ground::get<ball::velocity>(dd));
  siconos::remove(xball1, dd);
  print("{},{},{}\n", xball1.get(), xball2.get(), xball3.get());
  print("{}\n", ground::get<ball::velocity>(dd));

  static_assert(match::tag<attached_storage<ball, zz, some::scalar>, zz>);
  auto xball4 = add<ball>(dd);
  xball4.property(symbol<"x">{}) = 3.0;
  xball4.property("z"_s) = 4.0;
  xball2.property(symbol<"x">{}) = 1.0;
  xball2.property("z"_s) = 2.0;



  print("{},{},{},{}\n", xball2.property(symbol<"x">{}), xball2.property("z"_s), xball4.property(symbol<"x">{}), xball4.property("z"_s));

   for_each(dd, []<typename Key, typename Value>(Key k, Value v)
            {
              print("Key: [{}]\n", ground::type_name<Key>());
              print("Value: [{}]\n\n", ground::type_name<Value>());
            });
//   assert((ground::get<attached_storage<some::unbounded_collection<ball>, internal_index, some::indice>>(dd)[0][1] == 1));

}
