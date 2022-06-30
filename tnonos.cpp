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
      return siconos::add_item(step, data);
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
  template<typename Param>
  struct lagrangian
  {
    static constexpr auto dof = Param::dof;
    struct dynamical_system
    {
      struct mass_matrix
      {
        using type = types::matrix<dof, dof>;
      };

      struct q
      {
        using type = types::vector<dof>;
      };

      struct velocity
      {
        using type = types::vector<dof>;
      };

      struct fext
      {
        using type = types::vector<dof>;
      };

      using attributes = std::tuple<mass_matrix, q, velocity, fext>;

    };

    struct relation
    {
    };

    using items = std::tuple<dynamical_system, relation>;
  };

  template<typename Param>
  struct nonsmooth_law
  {

    struct newton_impact_friction
    {
      struct e
      {
        using type = types::scalar;
      };

      struct mu
      {
        using type = types::scalar;
      };

      using attributes = std::tuple<e, mu>;
    };

    struct newton_impact
    {
      struct e
      {
        using type = types::scalar;
      };

      using attributes = std::tuple<e>;
    };
  };


  struct one_step_integrator
  {
    template<typename Form, typename Param>
    struct moreau_jean
    {

      static constexpr auto theta = Param::theta;

      using formulation = Form;
      using system = typename formulation::dynamical_system;


      template<std::size_t N, typename ...As>
      struct keeper
      {
        static constexpr std::size_t value = N;
        using type = std::tuple<As...>;
      };

      using keep = keeper<2, typename system::q, typename system::velocity>;

    };
  };

  template<typename Param>
  struct time_discretization
  {
  };
}

struct moreau_jean_param
{
  static constexpr auto theta = 0.5;
};

using namespace siconos;

int main()
{
  using formulation = lagrangian<param>;
  using osi = one_step_integrator::moreau_jean<formulation, moreau_jean_param>;
  using dynamical_system = formulation::dynamical_system;
  using velocity_t = dynamical_system::velocity;
  using siconos::get;
  auto data = siconos::make_data<env, osi, dynamical_system>();

  print("---\n");
  for_each(
    [](auto& a)
    {
      print("{}", std::size(a));
    },
    data._collections);

  print("---\n");

  auto ds0 = add_item<dynamical_system>(0, data);
  auto& v0 = get<velocity_t>(ds0, 0, data);

  v0 = { 1., 2., 3.};

  auto& velocityp = get<dynamical_system::velocity>(ds0, 0, data);
  print("--->{}\n", velocityp);

  for_each([](auto& a) { print("{:d}\n", std::size(a)); }, data._collections);
  print("{0}\n", data._collections);

  for_each([](auto& a) { print("{0}\n", a); }, data._collections);

  auto ds1 = add_item<dynamical_system>(0, data);

  get<dynamical_system::q>(ds1, 0, data) = { 1.,1., 1.};
  auto ds2 = add_item<dynamical_system>(0, data);
  get<dynamical_system::q>(ds2, 0, data) = { 9.,9., 9.};
  print("---\n");
  for_each([](auto& a) { print("{0}\n", a); }, data._collections);

  remove_item(ds1, 0, data);
  print("---\n");
  for_each([](auto& a) { print("{0}\n", a); }, data._collections);

  print("{}", get<dynamical_system::mass_matrix>(ds0, 0, data));

  auto m = get<dynamical_system::mass_matrix>(ds0, 0, data);

  print("---\n");

  print("{}", m);

  print("---\n");

  get<dynamical_system::fext>(ds0, 0, data) = { 10., 0., 0.};
  print("{}\n", data._collections);

  data(get<dynamical_system::fext>)(ds0, 0) = { 2.,1.,0.};
  print("{}\n", data(get<dynamical_system::fext>)(ds0, 0));

}

