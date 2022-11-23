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

    template<typename V, typename E>
    using graph = SiconosGraph <V, E,
                                boost::no_property,
                                boost::no_property,
                                boost::no_property>;

    template<typename ...Ts>
    using tuple = std::tuple<Ts...>;

    template<typename S>
    using unbounded_collection = std::vector<S>;

    template<typename S, std::size_t N>
    using bounded_collection = std::array<S, N>;

    template<typename T, indice N>
    using vector = std::array<T, N>;

    template<typename T, indice N, indice M>
    using matrix = std::array<T, N*M>;

    template<typename T, indice N>
    using diagonal_matrix = std::array<T, N>;

    template<typename T>
    using item_ref = siconos::internal_handle<T, indice>;

    template<typename T>
    using default_storage = std::array<T, 1>;

  };
}
#endif

