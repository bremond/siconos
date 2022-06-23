#include "SiconosData.hpp"

#include <vector>
#include <array>
#include <fmt/core.h>
#include <fmt/ranges.h>
#include "SiconosGraph.hpp"

using fmt::print;

struct env
{
  struct scalar
  {
    using type = double;
  };

  struct indice
  {
    using type = size_t;
  };

  struct graph
  {
    using type = SiconosGraph < indice::type, indice::type,
                                boost::no_property,
                                boost::no_property,
                                boost::no_property >;
  };

  struct vdescriptor
  {
    using type = graph::type::VDescriptor;
  };

  template<typename T>
  struct collection
  {
    using type = std::vector<T>;
  };

  template<indice::type N>
  struct vector
  {
    using type = std::array<scalar::type, N>;
  };

  template<indice::type N, indice::type M>
  struct matrix
  {
    using type = std::array<scalar::type, N*M>;
  };

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

      struct mass_matrix
      {
        using type = typename env::matrix<dof,dof>::type;
      };

      struct q { using type = typename env::vector<dof>::type; };

      struct velocity { using type = typename env::vector<dof>::type; };

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
      struct e : scalar{ scalar::type value; } _e;
      struct mu : scalar{} _mu;

      using attributes = std::tuple<e, mu>;
    };

    struct newton_impact
    {
      struct e : scalar{};

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


      template<env::indice::type N, typename ...As>
      struct keeper
      {
        static constexpr env::indice::type value = N;
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
  static constexpr typename Env::scalar::type theta=0.5;
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

