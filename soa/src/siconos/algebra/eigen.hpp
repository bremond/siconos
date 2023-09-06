#pragma once

#include <Eigen/Dense>

#include "siconos/algebra/linear_algebra.hpp"

namespace siconos {

// concepts
namespace match {

template <typename T>
concept matrix = requires(T m) { m(0, 0); };

template <typename T>
concept vector = matrix<T> && T::ColsAtCompileTime == 1;

template <typename T>
concept diagonal_matrix = !matrix<T> && requires(T m) { m.diagonal()[0]; };

template <typename T>
concept any_matrix = (diagonal_matrix<T> || matrix<T>);

}  // namespace match

template <typename T, size_t M, size_t N>
using matrix = Eigen::Matrix<T, M, N>;  // column storage

template <typename T, size_t M>
using vector = Eigen::Vector<T, M>;  // column vector

template <typename T>
using matrix_view = Eigen::Map<T>;
static_assert(vector<int, 3>::ColsAtCompileTime == 1);

template <typename T, size_t M>
using diagonal_matrix = Eigen::DiagonalMatrix<T, M>;

template <typename T>
struct value_type {
  using type = decltype([]<bool flag = false>() {
    if constexpr (requires(T m) { typename T::value_type; }) {
      return typename T::value_type{};  // ok with gcc!
    }
    else {
      static_assert(flag, "no value_type");
    }
  }());
};

// template specialization ok with clang, fails with gcc:
//
// template <typename T, size_t M, size_t N>
// struct value_type<Eigen::Matrix<T, M, N>> {
//   using type = typename Eigen::template Matrix<T, M, N>::value_type;
// };

// template <typename T, size_t M>
// struct value_type<Eigen::DiagonalMatrix<T, M>> {
//   using type = typename Eigen::template DiagonalMatrix<T, M>::Scalar;
// };

template <typename A>
using trans_t = matrix<typename value_type<A>::type, A::ColsAtCompileTime,
                       A::RowsAtCompileTime>;

template <typename A, typename B>
using prod_t = matrix<typename value_type<B>::type, A::RowsAtCompileTime,
                      B::ColsAtCompileTime>;

static_assert(std::is_same_v<matrix<int, 1, 2>, trans_t<matrix<int, 2, 1>>>);
static_assert(std::is_same_v<prod_t<matrix<int, 2, 2>, matrix<int, 2, 1>>,
                             matrix<int, 2, 1>>);
static_assert(std::is_same_v<prod_t<matrix<int, 1, 2>, matrix<int, 2, 1>>,
                             matrix<int, 1, 1>>);
static_assert(
    std::is_same_v<prod_t<matrix<int, 1, 2>, trans_t<matrix<int, 1, 2>>>,
                   matrix<int, 1, 1>>);

}  // namespace siconos
