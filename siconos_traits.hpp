#ifndef SICONOS_TRAITS_HPP
#define SICONOS_TRAITS_HPP

#include "siconos_pattern.hpp"
namespace siconos
{
  namespace traits
  {

    template<typename E, match::attribute T>
    struct config
    {
      using type = decltype(
        []()
        {
          if constexpr (std::derived_from<T, some::scalar>)
          {
            return typename E::scalar{};
          }
          else if constexpr (std::derived_from<T, some::indice>)
          {
            return typename E::indice{};
          }
          else if constexpr (match::vdescriptor<T>)
          {
            return typename E::template vdescriptor<typename T::type>{};
          }
          else if constexpr (std::tuple_size_v<decltype(T::sizes)> == 1)
          {
            return typename E::template vector<std::get<0>(T::sizes)>{};
          }
          else if constexpr (std::tuple_size_v<decltype(T::sizes)> == 2)
          {
            return typename E::template matrix<
              std::get<0>(T::sizes),
              std::get<1>(T::sizes)>{};
          }
          else
          {
            // not found
            // cf https://stackoverflow.com/questions/38304847/constexpr-if-and-static-assert
            []<bool flag = false>()
              {
                static_assert(T{} == flag, "cannot translate type");
              }();
          }

        }());
    };
  }
}
#endif
