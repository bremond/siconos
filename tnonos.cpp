#include "SiconosData.hpp"

#include <vector>
#include <array>
#include <fmt/core.h>
#include <fmt/ranges.h>
#include "SiconosGraph.hpp"

using fmt::print;

struct env
{
  using scalar = double;
  using indice = size_t;

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
  static constexpr size_t dof = 3;
};

namespace siconos
{
  template<typename Param>
  struct lagrangian
  {
    using param = Param;

    static constexpr auto dof = param::dof;

    struct time_invariant_dynamical_system
    {
      struct mass_matrix
      {
        template<typename Env>
        using type = Env::template matrix<dof, dof>;
      } mass_matrix;

      struct q
      {
        template<typename Env>
        using type = Env::template vector<dof>;
      } q;

      struct velocity
      {
        template<typename Env>
        using type = Env::template vector<dof>;
      } velocity;
    };

    struct relation
    {
    };
  };

  template<typename Env, typename Param>
  struct nonsmooth_law
  {
    struct newton_impact_friction
    {
    };
    struct newton_impact
    {
      struct e
      {
        template<typename Env>
        using type = typename env::scalar;
      } e;
    };
  };


  template<typename Env>
  struct one_step_integrator
  {
    template<typename Formulation, typename Param>
    struct moreau_jean
    {
      using formulation = Formulation;

      static constexpr typename Formulation::env::indice memory_size = 2;
    }
  };

  template<typename Env>
  struct time_discretization
  {
  }
}

using namespace siconos;

int main()
{
  using formulation = lagrangian<param>;
  using osi = one_step_integrator<env>::moreau_jean<formulation>;

  auto data_ds = make_data<osi>();

  for_each([](auto& a) { print("{:d}\n", std::size(a)); }, data_ds.collections);

  add_item(data_ds);

  auto& velocity = get(ds().velocity, data_ds);
  velocity[0] = 1.0;
  velocity[1] = 2.0;

  print("{0}", data_ds.collections);

  for_each([](auto& a) { print("{0}\n", a); }, data_ds.collections);

  add_item(data_ds);
  add_item(data_ds);
  print("---\n");
  for_each([](auto& a) { print("{0}\n", a); }, data_ds.collections);

  remove_item(1, data_ds);
  print("---\n");
  for_each([](auto& a) { print("{0}\n", a); }, data_ds.collections);
}

