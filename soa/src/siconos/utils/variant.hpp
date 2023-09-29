#include <tuple>
#include <variant>
#include <siconos/storage/ground/ground.hpp>
namespace siconos::utils {

template <typename T>
struct unfold_variant_s;

template <typename... Ts>
struct unfold_variant_s<std::variant<Ts...>> {
  using types = storage::ground::tuple<Ts...>;

  static constexpr std::size_t number_of_types = sizeof...(Ts);

};

const constexpr auto& unfold_variant(auto& v)
{
  return unfold_variant_s<std::decay_t<decltype(v)>>{};
};

template<typename T, typename V>
V apply_maybe(T&& var, V&& default_result, auto&& fun)
{
  switch (var.index())
  {
  case 0: fun(std::get<0>(static_cast<T&&>(var)));
  default: return default_result;
  }
}
}
