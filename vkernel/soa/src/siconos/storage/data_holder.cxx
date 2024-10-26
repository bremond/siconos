module;

export module siconos.storage:data_holder;

import siconos.storage.ground;
import siconos.storage.pattern;
import siconos.storage.some;
import :info;
import :get;

export namespace siconos::storage {

using namespace pattern;

template <typename Handle>
struct default_interface {
  decltype(auto) self()
  {
    return static_cast<Handle*>(
        this);  // handle inherits from default_interface
  };

  auto env()
  {
    auto& data = self()->data();
    using info_t = get_info_t<decltype(data)>;
    return typename info_t::env{};
  }

  auto params() { return typename decltype(env())::params{}; }

  template <ground::string_literal S>
  constexpr auto env_param()
  {
    return ground::get<pattern::param<S>>(self()->params()).value;
  }
};

template <typename Struct>
struct data_holder : item<> {
  using attributes =
      gather<attribute<"instance", some::specific<pointer<Struct>>>>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) instance() { return attr<"instance">(*self()); };
  };
};

}  // namespace siconos::storage
