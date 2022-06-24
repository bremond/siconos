#include "SiconosData.hpp"

#include <vector>
#include <array>
#include <fmt/core.h>
#include <fmt/ranges.h>
#include "SiconosGraph.hpp"

using fmt::print;

template<typename X, typename T >
struct attribute
{
  using type = T;
  static constexpr T& get(const auto vd, const auto step, auto& data)
  { return siconos::get<X>(vd, step, data); };
};


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
  struct dof
  {
    static constexpr auto value = 3;
  };
};

namespace siconos
{

  template<typename Param>
  struct lagrangian
  {
    using param = Param;

    static constexpr auto dof = param::dof::value;

    template<typename Env>
    struct dynamical_system
    {
      using env = Env;

      struct mass_matrix :
        attribute<mass_matrix,
                  typename env::matrix<dof,dof>> {};

      struct q : attribute<q, typename env::vector<dof>>{};

      struct velocity : attribute<velocity, typename env::vector<dof>>{};

      using attributes = std::tuple<mass_matrix, q, velocity>;

    };

    template<typename Env>
    struct relation
    {
    };
  };

  template<typename Env, typename Param>
  struct nonsmooth_law
  {
    using env = Env;
    using scalar = env::scalar;

    struct newton_impact_friction
    {
      struct e : attribute<e, scalar>{} _e;
      struct mu : attribute<mu, scalar>{} _mu;

      using attributes = std::tuple<e, mu>;
    };

    struct newton_impact
    {
      struct e : attribute<e, scalar>{} _e;

      using attributes = std::tuple<e>;
    };
  };


  template<typename Env>
  struct one_step_integrator
  {
    template<typename Form, typename Param>
    struct moreau_jean
    {
      using env = Env;


      static constexpr auto theta = Param::theta;

      using formulation = Form;
      using system = typename formulation::dynamical_system<env>;


      template<env::indice N, typename ...As>
      struct keeper
      {
        static constexpr env::indice value = N;
        using type = std::tuple<As...>;
      };

      using keep = keeper<2, typename system::q, typename system::velocity>;

    };
  };

  template<typename Env>
  struct time_discretization
  {
  };
}

template<typename Env>
struct moreau_jean_param
{
  static constexpr typename Env::scalar theta=0.5;
};

using namespace siconos;

int main()
{
  using formulation = lagrangian<param>;
  using osi = one_step_integrator<env>::moreau_jean<formulation, moreau_jean_param<env>>;
  using dynamical_system = formulation::dynamical_system<env>;
  using velocity_t = dynamical_system::velocity;

  auto data = make_data<osi>();

  print("---\n");
  for_each([](auto& a) { print("{:d}\n", std::size(a)); }, data.collections);
  print("{}", data.collections);
  print("---\n");

  auto ds0 = add_item(0, data);
  auto& velocity = get<velocity_t>(ds0, 0, data);
  velocity[0] = 1.0;
  velocity[1] = 2.0;
  for_each([](auto& a) { print("{:d}\n", std::size(a)); }, data.collections);
  print("{0}\n", data.collections);

  for_each([](auto& a) { print("{0}\n", a); }, data.collections);

  auto ds1 = add_item(0, data);
  dynamical_system::q::get(ds1, 0, data) = std::array{ 1.,1., 1.};
  auto ds2 = add_item(0, data);
  dynamical_system::q::get(ds2, 0, data) = std::array{ 9.,9., 9.};
  print("---\n");
  for_each([](auto& a) { print("{0}\n", a); }, data.collections);

  remove_item(1, 1, data);
  print("---\n");
  for_each([](auto& a) { print("{0}\n", a); }, data.collections);

  print("{}", get<dynamical_system::mass_matrix>(ds0, 0, data));

  auto m = dynamical_system::mass_matrix::get(ds0, 0, data);

  print("---\n");

  print("{}", m);
}

