
#include <siconos/storage/pattern/pattern.hpp>

using namespace siconos::storage::pattern;

static_assert(
  std::is_same_v<decltype(cons(int{}, std::tuple<char, float>{})),
  std::tuple<int, char, float>>);

static_assert(
  std::is_same_v<decltype(append(gather<int>{}, gather<float>{})),
  gather<int, float>>);

static_assert(transform([]<typename T>(T) { return int{}; },
                        std::tuple<char, float, double>{}) ==
              std::tuple<int, int, int>{});



int main() {};
