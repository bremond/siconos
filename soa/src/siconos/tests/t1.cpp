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
  struct interface : siconos::default_interface<H> {};
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

static_assert(
    match::diagonal_matrix<
        decltype(storage::attr<"attr0">(storage::add<item0>(
            storage::make<standard_environment<int>, item0,
                          storage::with_properties<
                              storage::diagonal<item0, "attr0">>>())))>);

static_assert(
    match::matrix<decltype(storage::attr<"attr0">(storage::add<item0>(
        storage::make<standard_environment<int>,
                      wrap<some::unbounded_collection, item0>>())))>);

static_assert(
    match::diagonal_matrix<
        decltype(storage::attr<"attr0">(storage::add<item0>(
            storage::make<standard_environment<int>,
                          wrap<some::unbounded_collection, item0>,
                          storage::with_properties<
                              storage::diagonal<item0, "attr0">>>())))>);

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

static_assert(
    match::diagonal_matrix<
        decltype(storage::attr<"mass_matrix">(storage::add<ball>(
            storage::make<
                standard_environment<config::map<config::iparam<"dof", 3>>>,
                simulation, ball, relation, interaction,
                storage::with_properties<
                    storage::diagonal<ball, "mass_matrix">>>())))>);

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

namespace siconos::config {
using siconos::storage::pattern::with_name;
using disk = with_name<"disk", model::lagrangian_ds>;
using lcp =
    with_name<"lcp", simul::nonsmooth_problem<LinearComplementarityProblem>>;
using osnspb = with_name<"osnspb", simul::one_step_nonsmooth_problem<lcp>>;
using nslaw = with_name<"nslaw", model::newton_impact>;
using diskdisk_r = with_name<"diskdisk_r", model::diskdisk_r>;
using diskplan_r = with_name<"diskplan_r", model::diskplan_r>;
using interaction =
    with_name<"interaction",
              simul::interaction<nslaw, diskdisk_r, diskplan_r>>;
using osi =
    with_name<"osi",
              simul::one_step_integrator<disk, interaction>::moreau_jean>;
using td = with_name<"time_discretization", simul::time_discretization<>>;
using topo = with_name<"topology", simul::topology<disk, interaction>>;
using simulation =
    with_name<"simulation", simul::time_stepping<td, osi, osnspb, topo>>;

using params = map<iparam<"dof", 3>>;
}  // namespace siconos::config::disks

namespace pattern = siconos::storage::pattern;
int main()
{

  auto ddd = storage::make<env, with_name<"bbb", bbb>>();

  auto bob1 = storage::add<bbb>(ddd);

  bbb::attr::at(bob1).reset(new aaa);

  struct item1 : item<> {
    using attributes = gather<attribute<"one", some::scalar>>;
  };

  auto eee = storage::make<env, item1>();

  auto h1 = storage::add<item1>(eee);

  storage::attr<"one">(h1) = 1.0;

  // auto data = storage::make<
  //     standard_environment<config::params>, config::simulation,
  //     pattern::wrap<some::unbounded_collection, config::disk>,
  //     config::diskdisk_r, config::diskplan_r,
  //     pattern::wrap<some::unbounded_collection, config::interaction>,
  //     storage::with_properties<
  //         storage::attached<config::disk, storage::pattern::symbol<"shape">,
  //                           storage::some::item_ref<model::disk>>,
  //         storage::time_invariant<
  //             storage::pattern::attr_t<config::disk, "fext">>,
  //         storage::diagonal<config::disk, "mass_matrix">,
  //         storage::unbounded_diagonal<storage::pattern::attr_t<
  //             config::osi, "mass_matrix_assembled">>>>();

  auto data = storage::make<
      standard_environment<config::params>, config::simulation,
      pattern::wrap<some::unbounded_collection, config::disk>,
      config::diskdisk_r, config::diskplan_r,
    pattern::wrap<some::unbounded_collection, config::interaction>>();
      // storage::with_properties<
      //     storage::attached<config::disk, storage::pattern::symbol<"shape">,
      //                       storage::some::item_ref<model::disk>>,
      //     storage::time_invariant<
      //         storage::pattern::attr_t<config::disk, "fext">>,
      //     storage::diagonal<config::disk, "mass_matrix">,
      //     storage::unbounded_diagonal<storage::pattern::attr_t<
      //         config::osi, "mass_matrix_assembled">>>>();

  using idata_t = decltype(data);
  using disks_info_t = std::decay_t<decltype(ground::get<storage::info>(
                                              idata_t{}))>;

  using disks_items_t = typename disks_info_t::all_items_t;
  auto named_items =
    ground::filter(disks_items_t{},
                   ground::derive_from<pattern::any_symbol>);

  ground::for_each(named_items, [](auto item) { std::cout << pattern::item_name(item) << std::endl; });
  //ground::dump_keys(data, [](auto&& key) { std::cout << key << std::endl;});
}
