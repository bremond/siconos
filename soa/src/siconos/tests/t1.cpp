#include "siconos/siconos.hpp"
#include "siconos/utils/print.hpp"

using namespace siconos;
using env = siconos::standard_environment;

#include "siconos/storage/ground/ground.hpp"
#include "siconos/storage/some/some.hpp"
#include "siconos/storage/traits/traits.hpp"
#include "siconos/utils/environment.hpp"

namespace siconos {

namespace some = storage::some;
namespace traits = storage::traits;
namespace ground = storage::ground;
using namespace storage::pattern;

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

int main()
{
  using ball = model::lagrangian_ds;
  using lcp = simul::nonsmooth_problem<LinearComplementarityProblem>;
  using osnspb = simul::one_step_nonsmooth_problem<lcp>;
  using relation = model::lagrangian_r;
  using nslaw = model::newton_impact;
  using interaction = simul::interaction<nslaw, relation, 1>;
  using osi = simul::one_step_integrator<ball, interaction>::moreau_jean;
  using td = simul::time_discretization<>;
  using topo = simul::topology<ball, interaction>;
  using simulation = simul::time_stepping<td, osi, osnspb, topo>;

  struct zz {};
  //  {
  static_assert(must::contains<int, std::tuple<int, float, char>>);
  static_assert(!must::contains<double, std::tuple<int, float, char>>);
  static_assert(
      std::is_same_v<decltype(cons(int{}, std::tuple<char, float>{})),
                     std::tuple<int, char, float>>);

  static_assert(
      std::is_same_v<decltype(append(gather<int>{}, gather<float>{})),
                     gather<int, float>>);

  static_assert(transform([]<typename T>(T) { return int{}; },
                          std::tuple<char, float, double>{}) ==
                std::tuple<int, int, int>{});
  static_assert(std::is_same_v<decltype(all_items(nslaw{})),
                               gather<siconos::model::newton_impact>>);
  static_assert(
      std::is_same_v<decltype([]() {
                       return transform(
                           []<typename T>(T) { return typename T::type{}; },
                           filter<hold<decltype([]<typename T>(T) {
                             return match::item_ref<T>;
                           })>>(typename interaction::attributes{}));
                     }()),
                     gather<nslaw, relation>>);

  // ground::type_trace<decltype(all_items(interaction{}))>();
  static_assert(std::is_same_v<decltype(all_items(interaction{})),
                               std::tuple<interaction, nslaw, relation>>);

  // static_assert(must::contains<osnspb, decltype(all_items(simulation{}))>);

  static_assert(match::item<ball>);
  static_assert(match::attribute<nslaw::e>);

  static_assert(match::attribute_of<nslaw::e, nslaw>);
  static_assert(match::attribute_of<ball::velocity, ball>);
  static_assert(match::attribute_of<td::step, td>);

  // static_assert(std::is_same_v<td,
  // decltype(item_attribute<td::step>(all_items(simulation{})))>);

  static_assert(
      std::is_same_v<
          typename siconos::traits::config<env>::convert<
              typename siconos::simul::time_discretization<>::step>::type,
          typename env::indice>);

  static_assert(
      std::is_same_v<
          typename siconos::traits::config<env>::convert<some::scalar>::type,
          typename env::scalar>);

  static_assert(
      match::attached_storage<storage::attached<ball, zz, some::scalar>,
                              ball>);
  static_assert(
      match::attached_storage<storage::attached<ball, zz, some::scalar>,
                              wrap<some::unbounded_collection, ball>>);
  static_assert(match::attached_storage<
                storage::attached<ball, symbol<"Z">, some::scalar>,
                wrap<some::unbounded_collection, ball>>);
  static_assert(
      !match::attached_storage<storage::attached<ball, zz, some::scalar>,
                               nslaw>);

  static_assert(
      match::unbounded_storage<some::unbounded_collection<some::scalar>>);
  static_assert(
      match::bounded_storage<
          some::bounded_collection<some::scalar, some::indice_value<1>>>);
  static_assert(
      !match::unbounded_storage<
          some::bounded_collection<some::scalar, some::indice_value<1>>>);
  static_assert(
      !match::bounded_storage<some::unbounded_collection<some::scalar>>);
  static_assert(match::unbounded_storage<
                some::unbounded_diagonal_matrix<some::scalar>>);

  static_assert(match::wrap<wrap<some::unbounded_diagonal_matrix, ball>>);

  static_assert(traits::translatable<int, env>);

  std::cout << ground::type_name<traits::config<env>::convert<int>::type>();
  static_assert(std::is_same_v<traits::config<env>::convert<int>::type, int>);

  static_assert(
      traits::translatable<siconos::some::unbounded_collection<char[3]>,
                           env>);
  static_assert(filter<hold<decltype([]<typename T>(T) {
                  return std::floating_point<T>;
                })>>(std::tuple<char, int, double>{}) ==
                std::tuple<double>{});

  static_assert(ground::filter(std::tuple<char, int, double>{},
                               ground::compose(
                                   ground::trait<std::is_floating_point>,
                                   ground::typeid_)) == std::tuple<double>{});

  static_assert(std::is_same_v<decltype(ground::filter(
                                   std::tuple<nslaw::e, ball::velocity>{},
                                   ground::derive_from<some::scalar>)),
                               std::tuple<nslaw::e>>);

  static_assert(
      std::is_same_v<
          decltype(attributes(interaction{})),
          gather<some::indice_parameter<"dof">, some::indice_value<1>,
                 interaction::nonsmooth_law, interaction::relation,
                 interaction::h_matrix1, interaction::h_matrix2,
                 interaction::lambda, interaction::y, interaction::ydot>>);

  static_assert(
      std::is_same_v<decltype(all_attributes(interaction{})),
                     gather<some::indice_parameter<"dof">,
                            some::indice_value<1>, interaction::nonsmooth_law,
                            interaction::relation, interaction::h_matrix1,
                            interaction::h_matrix2, interaction::lambda,
                            interaction::y, interaction::ydot, nslaw::e>>);

  //}

  auto ddd = storage::make_storage<env, bbb>();

  auto bob1 = storage::add<bbb>(ddd);

  bbb::attr::at(bob1).reset(new aaa);

  using item1 = item<attribute<"one", some::scalar>>;
  auto eee = storage::make_storage<env, item1>();

  auto h1 = storage::add<item1>(eee);

  storage::attr<"one">(h1) = 1.0;
}
