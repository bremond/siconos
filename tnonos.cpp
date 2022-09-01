#include "siconos_storage.hpp"
#include "siconos.hpp"

#include <tuple>
#include <fmt/core.h>
#include <fmt/ranges.h>

using fmt::print;


struct env
{
  using scalar = double;
  using indice = std::size_t;

  using graph = SiconosGraph < indice, indice,
                               boost::no_property,
                               boost::no_property,
                               boost::no_property >;

  template<typename ...Ts>
  using tuple = std::tuple<Ts...>;

  template<typename T>
  using vdescriptor = tuple<graph::VDescriptor, T>;

  template<typename T>
  using collection = std::vector<T>;

  template<indice N>
  using vector = std::array<scalar, N>;

  template<indice N, indice M>
  using matrix = std::array<scalar, N*M>;

  template<typename T>
  using item_ref = std::tuple<T,indice>;
};

struct param
{
  static constexpr auto dof = 3;
};

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
  using simulation = time_stepping<td, osi, osnspb, ball, interaction>;
  using siconos::get;

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
                std::tuple<siconos::interaction<siconos::nonsmooth_law::newton_impact, siconos::lagrangian<param>::relation>, siconos::nonsmooth_law::newton_impact, siconos::lagrangian<param>::relation>>);

  static_assert(must::contains<osnspb, decltype(all_items(simulation{}))>);

  static_assert(match::item<ball>);
  static_assert(match::vertex_item<ball>);
  static_assert(match::attribute<nslaw::e>);

  static_assert(match::attribute_of<nslaw::e, nslaw>);
  static_assert(match::attribute_of<ball::velocity, ball>);
  static_assert(match::attribute_of<td::step, td>);

  static_assert(std::is_same_v<td, decltype(item_attribute<td::step>(all_items(simulation{})))>);

  static_assert(std::is_same_v<typename siconos::traits::config<env, typename siconos::time_discretization<param>::step>::type,
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

//  static_assert(flatten(std::tuple<siconos::time_stepping<siconos::time_discretization<param>, siconos::one_step_integrator<siconos::lagrangian<param>>::moreau_jean, siconos::one_step_nonsmooth_problem<siconos::lcp>, siconos::lagrangian<param>::dynamical_system, siconos::interaction<siconos::nonsmooth_law::newton_impact, siconos::lagrangian<param>::relation>>, siconos::time_discretization<param>, siconos::one_step_integrator<siconos::lagrangian<param>>::moreau_jean, siconos::one_step_nonsmooth_problem<siconos::lcp>, siconos::lagrangian<param>::dynamical_system, siconos::interaction<siconos::nonsmooth_law::newton_impact, siconos::lagrangian<param>::relation>, siconos::nonsmooth_law::newton_impact, siconos::lagrangian<param>::relation>{}) == false);

//  static_assert(all_items(typename machinery<simulation>::vertex{}) == false);
  auto data = siconos::make_storage<env, simulation>();
//   add<simulation>(data);

  print("---\n");
  for_each(
    [](auto& a)
    {
      print("->{}\n", a.size());
    },
    data._collections);

   print("---\n");

   auto ds0 = fix(add<ball>)(data);

   ball::q::get(ds0) = { 0., 0., 1.};

   auto& v0 = ball::velocity::get(ds0);

   v0 = { 1., 2., 3.};

   auto& velocityp = get<ball::velocity>(ds0.first, data);
   print("--->{}\n", velocityp);

   for_each([](auto& a) { print("{:d}\n", a.size()); }, data._collections);
//   print("{0}\n", data._collections);

//   for_each([](auto& a) { print("{0}\n", a); }, data._collections);

   auto ds1 = fix(add<ball>)(data);

   ball::q::get(ds1) = { 1.,1., 1.};

   auto ds2 = fix(add<ball>)(data);
   ball::q::get(ds2) = { 9.,9., 9.};
//   print("---\n");
//   for_each([](auto& a) { print("{0}\n", a); }, data._collections);

   fix_map(remove_vertex_item)(ds1, 0);
//   print("---\n");
//   for_each([](auto& a) { print("{0}\n", a); }, data._collections);

   print("{}", get<ball::mass_matrix>(ds0.first, data));

//   auto m = get<ball::mass_matrix>(ds0, data);

//   print("---\n");

//   print("{}", m);

//   print("---\n");

   ball::fext::get(ds0) = { 10., 0., 0.};
//   print("{}\n", data._collections);

   ball::fext::get(ds0) = { 2.,1.,0.};
//   print("{}\n", data(get<ball::fext>)(ds0));

   add<interaction>(data);
   add<nslaw>(data);
// //  add_attributes<nslaw>(data);

   auto& e = siconos::get_memory<nslaw::e>(data);
   e[0][0] = 0.7;

   print("{}\n", siconos::get_memory<nslaw::e>(data));

}
