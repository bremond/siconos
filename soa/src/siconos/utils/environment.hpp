#pragma once

#include <boost/container/vector.hpp>
#include <cstddef>
#include <array>
#include <vector>
#include <tuple>
#include <boost/container/static_vector.hpp>
#include "SolverOptions.h"

#include "siconos/utils/SiconosGraph.hpp" // modified for std::array
#include "siconos/utils/pattern.hpp"
#include "siconos/algebra/numerics.hpp"
#include "siconos/siconos.hpp"

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

    using boolean = bool;
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
    using matrix = std::array<std::array<T, M>, N>; // row major

    template<typename T>
    using unbounded_matrix = numerics::mat<T>;

    template<typename T>
    using unbounded_vector = numerics::vec<T>;

    template<typename T, indice N>
    using diagonal_matrix = std::array<T, N>;

    template<typename T>
    using unbounded_diagonal_matrix = numerics::diag_mat<T>;

    template<typename T>
    using item_ref = siconos::index<T, indice>;

    template<typename T>
    using default_storage = std::array<T, 1>;

    static constexpr auto white_color = boost::white_color;

    static constexpr auto gray_color = boost::gray_color;

    static constexpr auto black_color = boost::black_color;

    template<typename T>
    using vdescriptor = typename T::VDescriptor;

    using params = ground::map<ground::pair<param<"dof">, param_val<3>>>;

  };


}


