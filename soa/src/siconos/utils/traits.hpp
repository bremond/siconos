#pragma once

#include "siconos/algebra/eigen.hpp"
#include "siconos/utils/pattern.hpp"

namespace siconos {
namespace traits {
static auto translate = rec([]<typename E, typename T>(auto&& translate, E,
                                                       T) {
  if constexpr (!match::attribute<T>) {
    // return the type itself
    return T{};
  }
  else if constexpr (std::derived_from<T, some::boolean>) {
    return typename E::boolean{};
  }
  else if constexpr (std::derived_from<T, some::scalar>) {
    return typename E::scalar{};
  }
  else if constexpr (std::derived_from<T, some::indice>) {
    return typename E::indice{};
  }
  else if constexpr (std::derived_from<T, some::undefined_indice_parameter>) {
    return ground::get<param<T::name>>(typename E::params{});
  }
  else if constexpr (std::derived_from<T, some::undefined_indice_value>) {
    return T{};
  }

  else if constexpr (std::derived_from<
                         T, some::undefined_unbounded_collection>) {
    return typename E::template unbounded_collection<decltype(translate(
        E{}, typename T::type{}))>{};
  }
  else if constexpr (std::derived_from<T,
                                       some::undefined_bounded_collection>) {
    return typename E::template bounded_collection<
        decltype(translate(E{}, typename T::type{})),
        decltype(translate(E{}, nth_t<0, typename T::sizes>{}))::value>{};
  }
  else if constexpr (std::derived_from<T, some::undefined_vector>) {
    return typename E::template vector<
        decltype(translate(E{}, typename T::type{})),
        decltype(translate(E{}, nth_t<0, typename T::sizes>{}))::value>{};
  }
  else if constexpr (std::derived_from<T, some::undefined_diagonal_matrix> &&
                     std::derived_from<T, some::unbounded_storage>) {
    return typename E::template unbounded_diagonal_matrix<decltype(translate(
        E{}, typename T::type{}))>{};
  }
  else if constexpr (std::derived_from<T, some::undefined_diagonal_matrix>) {
    return typename E::template diagonal_matrix<
        decltype(translate(E{}, typename T::type{})),
        decltype(translate(E{}, nth_t<0, typename T::sizes>{}))::value>{};
  }
  else if constexpr (std::derived_from<T, some::undefined_matrix>) {
    using nrows = decltype(translate(E{}, nth_t<0, typename T::sizes>{}));
    using ncols = decltype(translate(E{}, nth_t<1, typename T::sizes>{}));
    return typename E::template matrix<decltype(translate(
                                           E{}, typename T::type{})),
                                       nrows::value, ncols::value>{};
  }
  else if constexpr (std::derived_from<T, some::undefined_unbounded_matrix>) {
    return typename E::template unbounded_matrix<decltype(translate(
        E{}, typename T::type{}))>{};
  }
  else if constexpr (std::derived_from<T, some::undefined_unbounded_vector>) {
    return typename E::template unbounded_vector<decltype(translate(
        E{}, typename T::type{}))>{};
  }
  else if constexpr (std::derived_from<T, some::undefined_graph>) {
    return typename E::template graph<
        decltype(translate(E{}, std::get<0>(typename T::types{}))),
        decltype(translate(E{}, std::get<1>(typename T::types{})))>{};
  }
  else if constexpr (std::derived_from<T, some::undefined_vdescriptor>) {
    return typename E::template vdescriptor<decltype(translate(
        E{}, typename T::type{}))>{};
  }
  else if constexpr (std::derived_from<T, some::item_ref<typename T::type>>) {
    return typename E::template item_ref<typename T::type>{};
  }
  else if constexpr (std::derived_from<T, some::given_type>) {
    return typename T::type{};
  }
  else {
    []<typename Attr = T, bool flag = false>()
    {
      static_assert(flag, "cannot translate attribute");
    }
    ();
  }
});

template <typename T, typename E>
concept translatable = requires(T t) { translate(E{}, t); };

template <typename E>
struct config {
  template <translatable<E> A>
  struct convert {
    using type = decltype(translate(E{}, A{}));
  };
};

template <typename T>
static constexpr decltype(auto) nrows(T)
{
  if constexpr (match::matrix<T> || match::diagonal_matrix<T>) {
    return T::RowsAtCompileTime; // specific to Eigen
  } else {
    []<typename Attr = T, bool flag = false>()
      {
        static_assert(flag, "no value type");
      }
    ();
  }
}

template <typename T>
static constexpr decltype(auto) ncols(T)
{
  if constexpr (match::matrix<T> || match::diagonal_matrix<T>) {
    return T::ColsAtCompileTime; // specific to Eigen
  } else {
    []<typename Attr = T, bool flag = false>()
      {
        static_assert(flag, "no value type");
      }
    ();
  }
}

template <typename T>
using value_type = typename T::value_type; // eigen ok.

}  // namespace traits

}  // namespace siconos
