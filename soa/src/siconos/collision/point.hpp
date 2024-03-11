#pragma once

#include "siconos/collision/collision_head.hpp"
#include "siconos/model/lagrangian_ds.hpp"

#include <concepts>
namespace siconos::collision {

// a point linked to an item (a dynamical system or a shape)
template <match::item Item>
struct point : item<> {
  using item_t = Item;
  using attributes = gather<
      // 3D coordinates, 2D => last value = 0.
      attribute<"coord", some::vector<some::scalar, some::indice_value<3>>>,
      attribute<"item", some::item_ref<item_t>>>;  // reverse link

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) coord() { return storage::attr<"coord">(*self()); };
    decltype(auto) item() { return storage::attr<"item">(*self()); };

    void update()
    {
      auto& data = self()->data();

      if constexpr (std::derived_from<item_t, model::lagrangian_ds>) {
        // one body / one point
        if constexpr (std::derived_from<item_t, model::lagrangian_ds>) {
          auto hbody = storage::handle(data, item());
          storage::attr<"coord">(*self())[0] = storage::attr<"q">(hbody)(0);
          storage::attr<"coord">(*self())[1] = storage::attr<"q">(hbody)(1);
          storage::attr<"coord">(*self())[2] = 0.; /* 2D */
        }
      }
    }
  };
};
}  // namespace siconos::collision
