#ifndef SICONOS_TRAITS_HPP
#define SICONOS_TRAITS_HPP

#include "siconos_pattern.hpp"

namespace siconos
{
  namespace traits
  {
    static auto translate =
      rec([]<typename E, match::attribute Attribute>
          (auto&& translate, E, Attribute)
          {
            if constexpr (std::derived_from<Attribute, some::scalar>)
            {
              return typename E::scalar{};
            }
            else if constexpr (std::derived_from<Attribute, some::indice>)
            {
              return typename E::indice{};
            }
            else if constexpr (std::derived_from<Attribute,
                               some::undefined_unbounded_collection>)
            {
              return typename E::template unbounded_collection<decltype(translate(E{}, typename Attribute::type{}))>{};
            }
            else if constexpr (std::derived_from<Attribute,
                               some::undefined_bounded_collection>)
            {
              return typename E::template bounded_collection<decltype(translate(E{}, typename Attribute::type{})),
                                                             std::get<0>(Attribute::sizes)>{};
            }
            else if constexpr (std::derived_from<Attribute,
                               some::undefined_vector>)
            {
              return typename E::template vector<decltype(translate(E{}, typename Attribute::type{})),
                                                 std::get<0>(Attribute::sizes)>{};
            }
            else if constexpr (std::derived_from<Attribute,
                               some::undefined_diagonal_matrix>)
            {
              return typename E::template diagonal_matrix<decltype(translate(E{}, typename Attribute::type{})),
                                                          std::get<0>(Attribute::sizes)>{};
            }
            else if constexpr (std::derived_from<Attribute,
                               some::undefined_matrix>)
            {
              return typename E::template matrix<
                decltype(translate(E{}, typename Attribute::type{})),
                std::get<0>(Attribute::sizes),
                std::get<1>(Attribute::sizes)>{};
            }
            else if constexpr (std::derived_from<Attribute, some::undefined_graph>)
            {
              return typename E::template graph<
                decltype(translate(E{}, std::get<0>(typename Attribute::types{}))),
                decltype(translate(E{}, std::get<1>(typename Attribute::types{})))>{};
            }
            else if constexpr (std::derived_from<Attribute,
                               some::item_ref<typename Attribute::type>>)
            {
              return typename E::template item_ref<typename Attribute::type>{};
            }
            else
            {
              // not found
              // cf https://stackoverflow.com/questions/38304847/constexpr-if-and-static-assert
              []<typename U = Attribute, bool flag = false>()
                {
                  static_assert(flag, "cannot translate this attribute");
                }();
            }
          });

    template<typename T, typename E>
    concept translatable = requires (T t) { translate(E{}, t); };

    template<typename E>
    struct config
    {
      template<translatable<E> A>
      struct convert
      {
        using type = decltype(translate(E{}, A{}));
      };
    };
  }
}
#endif
