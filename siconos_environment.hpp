#ifndef SICONOS_ENVIRONMENT_HPP
#define SICONOS_ENVIRONMENT_HPP

#include <boost/container/vector.hpp>
#include <cstddef>
#include <array>
#include <vector>
#include <tuple>
#include <boost/container/static_vector.hpp>
#include "SiconosGraph.hpp"
#include "siconos_pattern.hpp"

namespace siconos
{
  struct standard_environment
  {
    using scalar = double;
    using indice = std::size_t;

    using graph = SiconosGraph < indice, indice,
                                 boost::no_property,
                                 boost::no_property,
                                 boost::no_property >;

    template<typename ...Ts>
    using tuple = std::tuple<Ts...>;

    template<typename T>
    using vdescriptor = tuple<graph::VDescriptor, T>;

    template<typename T>
    using collection = std::vector<T>;

    template<indice N>
    using vector = std::array<scalar, N>;

    template<indice N, indice M>
    using matrix = std::array<scalar, N*M>;

    template<typename T>
    using item_ref = siconos::internal_handle<T, indice>;

  };

  struct fixed_environment
  {
    using scalar = double;
    using indice = std::size_t;

    using graph = SiconosGraph < indice, indice,
                                 boost::no_property,
                                 boost::no_property,
                                 boost::no_property >;

    template<typename ...Ts>
    using tuple = std::tuple<Ts...>;

    template<typename T>
    using vdescriptor = tuple<graph::VDescriptor, T>;

    template<typename T>
    using collection = boost::container::vector<T>;

    template<indice N>
    using vector = std::array<scalar, N>;

    template<indice N, indice M>
    using matrix = std::array<scalar, N*M>;

    template<typename T>
    using item_ref = siconos::internal_handle<T, indice>;
  };
}
#endif
