
#include <siconos/storage/pattern/pattern.hpp>

using namespace siconos::storage::pattern;
using namespace siconos::storage::some;

static_assert(std::is_same_v<decltype(cons(int{}, gather<char, float>{})),
                             gather<int, char, float>>);

static_assert(std::is_same_v<decltype(append(gather<int>{}, gather<float>{})),
                             gather<int, float>>);

static_assert(transform([]<typename T>(T) { return int{}; },
                        gather<char, float, double>{}) ==
              gather<int, int, int>{});

//static_assert(
//    attribute_name(siconos::storage::pattern::attribute<
//                   "hello", siconos::storage::some::scalar>{}) ==
//    "hello");

int main(){};
