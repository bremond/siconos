#include <concepts>

#include "siconos/config/config.hpp"
#include "siconos/siconos.hpp"
#include "siconos/storage/pattern/base.hpp"
#include "siconos/utils/print.hpp"

namespace siconos::config {
using params = map<assoc<param<"dof">, param_val<3>>>;
}

using namespace siconos;

using env = standard_environment<config::params>;

#include "siconos/model/lagrangian_r.hpp"
#include "siconos/model/model.hpp"
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
  struct attr : some::specific<pointer<aaa>>, siconos::access<attr> {};
  using attributes = gather<attr>;
  template <typename H>
  struct interface : siconos::default_interface<H> {
    constexpr int hello() { return 36; };

    auto methods()
    {
      return collect(method(symbol<"hello">{}, &interface<H>::hello));
    };
  };
};

static_assert(std::is_same_v<traits::config<standard_environment<int>>::
                                 convert<some::scalar>::type,
                             double>);
static_assert(std::is_same_v<traits::config<standard_environment<int>>::
                                 convert<pointer<float>>::type,
                             pointer<float>>);

static_assert(
    std::is_same_v<traits::config<standard_environment<int>>::convert<
                       some::specific<pointer<float>>>::type,
                   pointer<float>>);

static_assert(
    match::abstract_matrix<
        attribute<"attr0", some::matrix<some::scalar, some::indice_value<1>,
                                        some::indice_value<1>>>>);

struct item0 : item<> {
  using attributes = gather<
      attribute<"attr0", some::matrix<some::scalar, some::indice_value<1>,
                                      some::indice_value<1>>>>;
};

static_assert(match::diagonal_matrix<
              decltype(storage::attr<"attr0">(storage::add<item0>(
                  storage::make<standard_environment<int>, item0,
                                storage::with_properties<storage::diagonal<
                                    attr_t<item0, "attr0">>>>())))>);

static_assert(
    match::matrix<decltype(storage::attr<"attr0">(storage::add<item0>(
        storage::make<standard_environment<int>,
                      wrap<some::unbounded_collection, item0>>())))>);

static_assert(match::diagonal_matrix<
              decltype(storage::attr<"attr0">(storage::add<item0>(
                  storage::make<standard_environment<int>,
                                wrap<some::unbounded_collection, item0>,
                                storage::with_properties<storage::diagonal<
                                    attr_t<item0, "attr0">>>>())))>);

static_assert(
    std::is_same_v<
        std::decay_t<decltype((storage::add<model::newton_impact>(
                                   storage::make<standard_environment<int>,
                                                 model::newton_impact>()))
                                  .e())>,
        typename standard_environment<int>::scalar>);

static_assert(
    std::is_same_v<
        std::decay_t<
            decltype(storage::add<collision::disk_shape>(
                         storage::make<
                             standard_environment<int>, collision::disk_shape,
                             storage::with_properties<storage::bind<
                                 collision::disk_shape, "disk_shape">>>())
                         .radius())>,
        typename standard_environment<int>::scalar>);

static_assert(storage::has_property<collision::disk_shape,
                                    storage::property::bind>(
    storage::make<standard_environment<int>, collision::disk_shape,
                  storage::with_properties<
                      storage::bind<collision::disk_shape, "disk_shape">>>()));

}  // namespace siconos

using namespace boost::hana::literals;

using ball = model::lagrangian_ds;
using lcp = simul::nonsmooth_problem<LinearComplementarityProblem>;
using osnspb = simul::one_step_nonsmooth_problem<lcp>;
using nslaw = model::newton_impact;
using relation = model::lagrangian_r<nslaw::size>;
using interaction = simul::interaction<nslaw, relation>;
using osi = simul::one_step_integrator<ball, interaction>::moreau_jean;
using td = simul::time_discretization<>;
using topo = simul::topology<ball, interaction>;
using simulation = simul::time_stepping<td, osi, osnspb, topo>;
using inter_manager = simul::interaction_manager<nslaw>;
using disk_shape = collision::disk_shape;

// static_assert(
//   storage::attr<"nslaws">(
//     storage::add<inter_manager>(
//         storage::make <
//         standard_environment<config::map<config::iparam<"ncgroups", 1>>>,
//         inter_manager, nslaw>())
//         .insert_nonsmooth_law(
//             storage::add<nslaw>(
//                 storage::make <
//                 standard_environment<config::map<config::iparam<"ncgroups", 1>>>,
//                 inter_manager, nslaw>()),
//             0, 0))(0, 0) == storage::add<nslaw>(
//                 storage::make <
//                 standard_environment<config::map<config::iparam<"ncgroups", 1>>>,
//                 inter_manager, nslaw>()));

static_assert(
    match::diagonal_matrix<
        decltype(storage::attr<"mass_matrix">(storage::add<ball>(
            storage::make<
                standard_environment<config::map<config::iparam<"dof", 3>>>,
                simulation, ball, relation, interaction,
                storage::with_properties<
                    storage::diagonal<attr_t<ball, "mass_matrix">>>>())))>);

static_assert(
    match::index<std::decay_t<decltype(storage::prop<
                                       "shape">(storage::add<ball>(
        storage::make<
            standard_environment<config::map<config::iparam<"dof", 3>>>, ball,
            storage::with_properties<storage::attached<
                ball, symbol<"shape">, some::item_ref<disk_shape>>>>())))>>);

template <typename T>
struct is_polymorhic : std::integral_constant<bool, []() {
  return match::polymorphic_type<T>;
}()> {};

static_assert(match::relation1<storage::handle<
                  relation, int, decltype(storage::make<env, relation>())>>);

//  {
static_assert(std::is_same_v<decltype(all_items(nslaw{})),
                             gather<siconos::model::newton_impact>>);

static_assert(
    std::derived_from<std::decay_t<decltype(std::get<0>(ground::filter(
                          typename interaction::attributes{},
                          ground::compose(ground::trait<is_polymorhic>,
                                          ground::typeid_))))>,
                      some::polymorphic_attribute<some::item_ref<relation>>>);

// ground::type_trace<decltype(all_items(interaction{}))>();
static_assert(std::is_same_v<decltype(all_items(interaction{})),
                             std::tuple<interaction, nslaw, relation>>);

// static_assert(must::contains<osnspb,
// decltype(all_items(simulation{}))>);

static_assert(match::item<ball>);
static_assert(match::attribute<attr_t<nslaw, "e">>);

static_assert(match::attribute_of<attr_t<nslaw, "e">, nslaw>);

static_assert(match::attribute_of<attr_t<ball, "velocity">, ball>);
static_assert(match::attribute_of<attr_t<td, "step">, td>);

// static_assert(std::is_same_v<td,
// decltype(item_attribute<td::step>(all_items(simulation{})))>);

static_assert(
    std::is_same_v<typename siconos::traits::config<env>::convert<attr_t<
                       siconos::simul::time_discretization<>, "step">>::type,
                   typename env::indice>);

static_assert(
    std::is_same_v<
        typename siconos::traits::config<env>::convert<some::scalar>::type,
        typename env::scalar>);

struct zz {};

static_assert(
    match::attached_storage<storage::attached<ball, zz, some::scalar>, ball>);
static_assert(
    match::attached_storage<storage::attached<ball, zz, some::scalar>,
                            wrap<some::unbounded_collection, ball>>);
static_assert(match::attached_storage<
              storage::attached<ball, symbol<"Z">, some::scalar>,
              wrap<some::unbounded_collection, ball>>);
static_assert(!match::attached_storage<
              storage::attached<ball, zz, some::scalar>, nslaw>);

static_assert(
    match::unbounded_storage<some::unbounded_collection<some::scalar>>);
static_assert(match::bounded_storage<
              some::bounded_collection<some::scalar, some::indice_value<1>>>);
static_assert(!match::unbounded_storage<
              some::bounded_collection<some::scalar, some::indice_value<1>>>);
static_assert(
    !match::bounded_storage<some::unbounded_collection<some::scalar>>);
static_assert(
    match::unbounded_storage<some::unbounded_diagonal_matrix<some::scalar>>);

static_assert(match::wrap<wrap<some::unbounded_diagonal_matrix, ball>>);

static_assert(match::bounded_storage<
              wrap<some::bounded_collection, ball,
                   some::indice_value<1>>::template wrapper<some::scalar>>);

static_assert(traits::translatable<int, env>);

static_assert(std::is_same_v<traits::config<env>::convert<int>::type, int>);

static_assert(
    traits::translatable<siconos::some::unbounded_collection<char[3]>, env>);
static_assert(filter<hold<decltype([]<typename T>(T) {
                return std::floating_point<T>;
              })>>(std::tuple<char, int, double>{}) == std::tuple<double>{});

static_assert(
    ground::filter(std::tuple<char, int, double>{},
                   ground::compose(ground::trait<std::is_floating_point>,
                                   ground::typeid_)) == std::tuple<double>{});

static_assert(std::is_same_v<
              decltype(ground::filter(
                  std::tuple<attr_t<nslaw, "e">, attr_t<ball, "velocity">>{},
                  ground::derive_from<some::scalar>)),
              std::tuple<attr_t<nslaw, "e">>>);

static_assert(
    std::is_same_v<
        decltype(attributes(interaction{})),
        gather<attr_t<interaction, "relation">, attr_t<interaction, "nslaw">,
               attr_t<interaction, "h_matrix1">,
               attr_t<interaction, "h_matrix2">,
               attr_t<interaction, "lambda">, attr_t<interaction, "y">,
               attr_t<interaction, "ydot">>>);

static_assert(
    std::is_same_v<
        decltype(all_attributes(interaction{})),
        gather<attr_t<interaction, "relation">, attr_t<interaction, "nslaw">,
               attr_t<interaction, "h_matrix1">,
               attr_t<interaction, "h_matrix2">,
               attr_t<interaction, "lambda">, attr_t<interaction, "y">,
               attr_t<interaction, "ydot">, attr_t<nslaw, "e">,
               attr_t<relation, "b">, attr_t<relation, "h_matrix">>>);

//}

int main() {}
