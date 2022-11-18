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
#include "siconos.hpp"

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

    template<typename S>
    using unbounded_collection = std::vector<S>;

    template<typename S, std::size_t N>
    using bounded_collection = std::array<S, N>;

    template<indice N>
    using vector = std::array<scalar, N>;

    template<indice N, indice M>
    using matrix = std::array<scalar, N*M>;

    template<indice N>
    using diagonal_matrix = std::array<scalar, N>;

    template<typename T>
    using item_ref = siconos::internal_handle<T, indice>;

    template<typename T>
    using default_storage = std::array<T, 1>;

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
    using unbounded_collection = boost::container::vector<T>;

    template<typename T, std::size_t N>
    using bounded_collection = std::array<T, N>;

    template<typename T>
    using default_storage = std::array<T, 1>;

    template<indice N>
    using vector = std::array<scalar, N>;

    template<indice N, indice M>
    using matrix = std::array<scalar, N*M>;

    template<typename T>
    using item_ref = siconos::internal_handle<T, indice>;
  };
}
#endif
