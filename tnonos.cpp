#include "SiconosData.hpp"

#include <vector>
#include <array>
#include <fmt/core.h>
#include <fmt/ranges.h>
#include "SiconosGraph.hpp"

using fmt::print;

template<typename T>
struct member
{
  using type = T;
  T value;
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

      struct mass_matrix : member<typename env::matrix<dof,dof>> {}
        _mass_matrix;

      struct q : member<typename env::vector<dof>>{} _q;

      struct velocity : member<typename env::vector<dof>>{} _velocity;

      using attributes =
        std::tuple<mass_matrix, q, velocity>;

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
      struct e : member<scalar>{} _e;
      struct mu : member<scalar>{} _mu;

      using attributes = std::tuple<e, mu>;
    };

    struct newton_impact
    {
      struct e : member<scalar>{} _e;

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

  using velocity_t = formulation::dynamical_system<env>::velocity;

  auto data_ds = make_data<osi>();

  for_each([](auto& a) { print("{:d}\n", std::size(a)); }, data_ds.collections);

  add_item(0, data_ds);

  auto& velocity = get<velocity_t>(0, data_ds);
  velocity[0] = 1.0;
  velocity[1] = 2.0;

  print("{0}", data_ds.collections);

  for_each([](auto& a) { print("{0}\n", a); }, data_ds.collections);

  add_item(0, data_ds);
  add_item(0, data_ds);
  print("---\n");
  for_each([](auto& a) { print("{0}\n", a); }, data_ds.collections);

  remove_item(0, 1, data_ds);
  print("---\n");
  for_each([](auto& a) { print("{0}\n", a); }, data_ds.collections);
}

