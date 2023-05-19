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


#include "siconos/utils/environment.hpp"
#include "siconos/utils/some.hpp"

namespace siconos {
struct aaa {
  int v = 1;
};
struct bbb : item<> {
  struct attr : some::specific<pointer<aaa>>,
                siconos::access<attr>,
                some::attribute<> {};
  using attributes = gather<attr>;
  template <typename H>
  struct interface : siconos::default_interface<H> {
  };
};

static_assert(
    std::is_same_v<
        traits::config<standard_environment>::convert<some::scalar>::type,
        double>);
static_assert(
    std::is_same_v<
        traits::config<standard_environment>::convert<pointer<float>>::type,
        pointer<float>>);

static_assert(std::is_same_v<traits::config<standard_environment>::convert<
                                 some::specific<pointer<float>>>::type,
                             pointer<float>>);

}  // namespace siconos



using namespace boost::hana::literals;

int main() {
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
                gather<some::indice_parameter<"dof">, some::indice_value<1>,interaction::nonsmooth_law, interaction::relation, interaction::h_matrix, interaction::lambda, interaction::y, interaction::ydot>>);

  static_assert(std::is_same_v<decltype(all_attributes(interaction{})),
                gather<some::indice_parameter<"dof">, some::indice_value<1>,interaction::nonsmooth_law, interaction::relation, interaction::h_matrix,
                interaction::lambda, interaction::y, interaction::ydot, nslaw::e>>);

  //}

   auto ddd = make_storage<env, bbb>();

   auto bob1 = add<bbb>(ddd);

   bbb::attr::at(bob1).reset(new aaa);
}
