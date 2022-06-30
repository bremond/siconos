#include "SiconosData.hpp"

#include <vector>
#include <array>
#include <fmt/core.h>
#include <fmt/ranges.h>
#include "SiconosGraph.hpp"

using fmt::print;

namespace siconos
{
  template<typename X>
  struct item
  {
    using type = X;

    static constexpr auto add(const auto step, auto& data)
    {
      return siconos::add_vertex_item(step, data);
    };
  };


  template<typename T>
  struct attribute
  {
    using type = T;

    static constexpr auto get = [](const auto vd, const auto step, auto& data) -> typename T::type&
    {
//      return siconos::get<T>(vd, step, data);
    };
  };

}

struct env
{
  using scalar = double;
  using indice = std::size_t;

  using graph = SiconosGraph < indice, indice,
                               boost::no_property,
                               boost::no_property,
                               boost::no_property >;

  using vdescriptor = graph::VDescriptor;

  template<typename T>
  using collection = std::vector<T>;

  template<indice N>
  using vector = std::array<scalar, N>;

  template<indice N, indice M>
  using matrix = std::array<scalar, N*M>;

};

struct param
{
  static constexpr auto dof = 3;
};

namespace siconos
{
  template<typename ...Args>
  using tuple = std::tuple<Args...>;

  template<typename Param>
  struct lagrangian
  {
    static constexpr auto dof = Param::dof;
    struct dynamical_system
    {
      struct mass_matrix
      {
        using type = any::matrix<dof, dof>;
      };

      struct q
      {
        using type = any::vector<dof>;
      };

      struct velocity
      {
        using type = any::vector<dof>;
      };

      struct fext
      {
        using type = any::vector<dof>;
      };

      using attributes = tuple<mass_matrix, q, velocity, fext>;

    };

    struct relation
    {
    };

    using vertex_items = tuple<dynamical_system>;
    using edge_items = tuple<relation>;
  };

  struct nonsmooth_law
  {

    struct newton_impact_friction
    {
      struct e
      {
        using type = any::scalar;
      };

      struct mu
      {
        using type = any::scalar;
      };

      using attributes = tuple<e, mu>;
    };

    struct newton_impact
    {
      struct e
      {
        using type = any::scalar;
      };

      using attributes = tuple<e>;
    };

    using items = tuple<newton_impact_friction, newton_impact>;
  };

  template<typename Form>
  struct one_step_integrator
  {
    struct moreau_jean
    {
      using formulation = Form;
      using system = typename formulation::dynamical_system;

      template<std::size_t N, typename ...As>
      struct keeper
      {
        static constexpr std::size_t value = N;
        using type = tuple<As...>;
      };

      using keep = keeper<2, typename system::q, typename system::velocity>;

      struct theta
      {
        using type = any::scalar;
      };

      using attributes = tuple<theta>;
    };

    using items = tuple<moreau_jean>;
  };

  template<typename Param>
  struct time_discretization
  {
  };
}

using namespace siconos;

int main()
{
  using formulation = lagrangian<param>;
  using ball = formulation::dynamical_system;
  using osi = one_step_integrator<formulation>::moreau_jean;
  using nslaw = nonsmooth_law::newton_impact;
  using siconos::get;
  auto data = siconos::make_data<env, osi, ball, nslaw>();

  print("---\n");
  for_each(
    [](auto& a)
    {
      print("{}", std::size(a));
    },
    data._collections);

  print("---\n");

  auto ds0 = add_vertex_item<ball>(0, data);
  auto& v0 = get<ball::velocity>(ds0, 0, data);

  v0 = { 1., 2., 3.};

  auto& velocityp = get<ball::velocity>(ds0, 0, data);
  print("--->{}\n", velocityp);

  for_each([](auto& a) { print("{:d}\n", std::size(a)); }, data._collections);
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

  get<ball::fext>(ds0, 0, data) = { 10., 0., 0.};
  print("{}\n", data._collections);

  data(get<ball::fext>)(ds0, 0) = { 2.,1.,0.};
  print("{}\n", data(get<ball::fext>)(ds0, 0));

  add_item<nslaw>(0, data);

  auto& e = siconos::get_array<nslaw::e>(data);
  e[0][0] = 0.7;

  print("{}\n", siconos::get_array<nslaw::e>(data));
}

