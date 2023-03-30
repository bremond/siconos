#ifndef SICONOS_ENVIRONMENT_HPP
#define SICONOS_ENVIRONMENT_HPP

#include <boost/container/vector.hpp>
#include <cstddef>
#include <array>
#include <vector>
#include <tuple>
#include <boost/container/static_vector.hpp>
#include "SiconosGraph.hpp"
#include "SolverOptions.h"
#include "siconos_pattern.hpp"
#include "siconos.hpp"

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/functional.hpp>
#include <boost/numeric/ublas/io.hpp>


namespace siconos
{
  namespace ublas = boost::numeric::ublas;
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
    using vector = ublas::fixed_vector<T, N>;

    template<typename T, indice N, indice M>
    using matrix = ublas::fixed_matrix<T, M, N>;

    template<typename T>
    using unbounded_matrix = ublas::compressed_matrix<T>;

    template<typename T, indice N>
    using diagonal_matrix = ublas::diagonal_matrix<T, ublas::row_major, ublas::bounded_array<T, N>>;

    template<typename T>
    using unbounded_diagonal_matrix = ublas::diagonal_matrix<T, ublas::row_major, ublas::unbounded_array<T>>;

    template<typename T>
    using item_ref = siconos::index<T, indice>;

    template<typename T>
    using default_storage = std::array<T, 1>;

    static constexpr auto white_color = boost::white_color;

    static constexpr auto gray_color = boost::gray_color;

    static constexpr auto black_color = boost::black_color;

    template<typename T>
    using vdescriptor = typename T::VDescriptor;

  };
}
#endif

