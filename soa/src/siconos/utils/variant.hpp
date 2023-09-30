#include <array>
#include <siconos/storage/ground/ground.hpp>
#include <tuple>
#include <variant>
namespace siconos::utils {

template <typename Variant, typename V, typename F>
  requires requires {
    {
      std::variant_size_v<Variant> + 0
    } -> std::same_as<std::size_t>;
  }
V apply_if_valid(Variant& var, V&& default_value, F&& fun)
{
  return siconos::storage::ground::call_with_index<
      std::variant_size_v<Variant>, V>(
      var.index(), static_cast<V&&>(default_value),
      [&var, &fun](auto i) { return fun(std::get<i>(var)); });
}

}  // namespace siconos::utils
