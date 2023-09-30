#include <siconos/storage/ground/ground.hpp>

using namespace siconos::storage;


struct A {
  using a_t = void;
};
struct B : A {};
struct C : B {};
struct D {};
template <typename T>
struct E : T {};

template <typename T>
concept match_a_t = requires { typename T::a_t; };

static_assert(std::is_same_v<decltype(ground::filter(std::tuple<A, B, C, D, E<B>>{},
                                                     ground::derive_from<A>)),
              std::tuple<A, B, C, E<B>>>);

static_assert(match_a_t<A>);
static_assert(match_a_t<B>);

static_assert(
    std::is_same_v<decltype(ground::filter(
                              std::tuple<A, B, C, D, E<B>>{},
                       ground::is_a_model<[]<typename T>() consteval {
                         return match_a_t<T>;
                       }>)),
    std::tuple<A, B, C, E<B>>>);

int main() {}
