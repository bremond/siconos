#pragma once

#include "siconos/utils/SiconosGraph.hpp" // modified for std::array
#include "siconos/utils/pattern.hpp"
#include "siconos/storage/storage.hpp"
#include "siconos/algebra/linear_algebra.hpp"
#include "siconos/algebra/eigen.hpp"
#include "siconos/algebra/numerics.hpp"


namespace siconos
{
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

    template<typename T, indice M>
    using vector = siconos::vector<T, M>;

    template<typename T, indice M, indice N>
    using matrix = siconos::matrix<T, M, N>;

    template<typename T, indice N>
    using diagonal_matrix = siconos::diagonal_matrix<T, N>;

    template<typename T>
    using unbounded_matrix = numerics::mat<T>;

    template<typename T>
    using unbounded_vector = numerics::vec<T>;


    template<typename T>
    using unbounded_diagonal_matrix = numerics::diag_mat<T>;

    template<typename T>
    using item_ref = siconos::storage::index<T, indice>;

    template<typename T>
    using default_storage = std::array<T, 1>;

    static constexpr auto white_color = boost::white_color;

    static constexpr auto gray_color = boost::gray_color;

    static constexpr auto black_color = boost::black_color;

    template<typename T>
    using vdescriptor = typename T::VDescriptor;

    // should go in a config struct
    using params = ground::map<ground::pair<param<"dof">, param_val<3>>>;

  };
}


