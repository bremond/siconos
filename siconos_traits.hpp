#ifndef SICONOS_TRAITS_HPP
#define SICONOS_TRAITS_HPP

#include "siconos_pattern.hpp"
namespace siconos
{
  namespace traits
  {
    struct null_type {};

    template<some::type T>
    struct missing_conversion_for
    {
      static constexpr auto failure =
        []<bool flag = false>
      {
        static_assert(flag);
      };
    };

    template<typename E, some::type T>
    struct config
    {
      using type = missing_conversion_for<T>;
    };

    template<typename E>
    struct config<E, some::scalar>
    {
      using type = typename E::scalar;
    };

    template<typename E>
    struct config<E, some::indice>
    {
      using type = typename E::indice;
    };

    template<typename E, typename T>
    struct config<E, some::vdescriptor<T>>
    {
      using type = typename E::vdescriptor<T>;
    };

    template<std::size_t N, typename E>
    struct config<E, some::vector<N>>
    {
      using type = typename E::template vector<N>;
    };

    template<std::size_t N, std::size_t M, typename E>
    struct config<E, some::matrix<N, M>>
    {
      using type = typename E::template matrix<N, M>;
    };
  };
}
#endif
