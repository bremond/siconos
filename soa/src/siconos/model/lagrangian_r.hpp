#pragma once

#include "siconos/model/model.hpp"

#include "siconos/algebra/algebra.hpp"
#include "siconos/model/lagrangian_ds.hpp"

namespace siconos::storage::pattern::match {
template <typename T>
concept linear_relation = match::handle<T, model::linear>;

template <typename T>
concept relation1 = match::handle<T, model::relation1>;

template <typename T>
concept relation2 = match::handle<T, model::relation2>;

}  // namespace siconos::storage::pattern::match

namespace siconos::model {

template <auto NSLSize>
struct lagrangian_r : item<>,
                      linear,
                      relation1,
                      relation2,
                      any_lagrangian_relation {
  using nslaw_size = some::param_val<NSLSize>;
  using dof = some::indice_parameter<"dof">;

  struct h_matrix : some::matrix<some::scalar, nslaw_size, dof>,
                    access<h_matrix> {};

  struct b : some::scalar, access<b> {};

  using attributes = gather<b, h_matrix>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) h_matrix() { return Handle::type::h_matrix::at(*self()); }

    decltype(auto) b() { return Handle::type::b::at(*self()); }

    decltype(auto) compute_jachq(auto step, auto& q1, auto& q2,
                                 auto& h_matrix1, auto& h_matrix2)
    {
      h_matrix1 = h_matrix();
      h_matrix2 = -h_matrix();
    }
    decltype(auto) compute_jachq(auto step, auto& q, auto& h_matrix1)
    {
      h_matrix1 << -h_matrix();
    }
  };
};

struct disk : item<> {
  using attributes = gather<attribute<"radius", some::scalar>>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) radius() { return attr<"radius">(*self()); };

    decltype(auto) distance(match::vector auto& q)
    {
      auto& x = q(0);
      auto& y = q(1);
    }
  };
};

struct plan : item<> {
  // a*x + b*y + c = 0
  using attributes = gather<
      attribute<"coeffs", some::vector<some::scalar, some::param_val<4>>>>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) a() { return attr<"coeffs">(*self())(0); };
    decltype(auto) b() { return attr<"coeffs">(*self())(1); };
    decltype(auto) c() { return attr<"coeffs">(*self())(2); };

    // 1/sqrt(a^2+b^2)
    decltype(auto) invsqrta2pb2() { return attr<"coeffs">(*self())(3); };

    decltype(auto) distance(match::vector auto& q)
    {
      auto& x = q(0);
      auto& y = q(1);

      return fabs(a() * x + b() * y + c()) * invsqrta2pb2();
    }
  };
};

  struct diskplan_r : item<>, relation1, any_lagrangian_relation {
  using attributes = gather<attribute<"disk", some::item_ref<disk>>,
                            attribute<"plan", some::item_ref<plan>>>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) r()
    {
      return handle(self()->data(), attr<"disk">(*self())).radius();
    };
    decltype(auto) plan()
    {
      return handle(self()->data(), attr<"plan">(*self()));
    };

    decltype(auto) compute_h(match::vector auto& q)
    {
      return plan().distance(q);
    }

    decltype(auto) compute_jachq(auto step, auto& q, auto& h_matrix1)
    {
      auto& x = q(0);
      auto& y = q(1);
      auto& a = plan().a();
      auto& b = plan().b();
      auto& c = plan().c();
      auto& invsqrta2pb2 = plan().invsqrta2pb2();

      auto d1 = a * x + b * y + c;
      auto sd1 = copysign(1, d1);

      h_matrix1(0, 0) = a * sd1 * invsqrta2pb2;
      h_matrix1(1, 0) = -b * sd1 * invsqrta2pb2;
      h_matrix1(0, 1) = b * sd1 * invsqrta2pb2;
      h_matrix1(1, 1) = a * sd1 * invsqrta2pb2;
      h_matrix1(0, 2) = 0;
      h_matrix1(1, 2) = -r();
    }
  };
};

  struct diskdisk_r : item<>, relation2, any_lagrangian_relation {
  using attributes = gather<attribute<"disk1", some::item_ref<disk>>,
                            attribute<"disk2", some::item_ref<disk>>>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) r1()
    {
      return handle(self()->data(), attr<"disk1">(*self())).radius();
    };
    decltype(auto) r2()
    {
      return handle(self()->data(), attr<"disk2">(*self())).radius();
    };

    decltype(auto) compute_h(match::vector auto& q1, match::vector auto& q2)
    {
      auto& dx = q2[0] - q1[0];
      auto& dy = q2[1] - q2[1];

      return algebra::hypot(dx, dy) - (r1() + r2());
    }

    decltype(auto) compute_jachq(auto step, auto& q1, auto& q2,
                                 auto& h_matrix1, auto& h_matrix2)
    {
      auto x1 = q1(0);
      auto y1 = q1(1);
      auto x2 = q2(0);
      auto y2 = q2(1);

      auto dx = x2 - x1;
      auto dy = y2 - y1;

      auto d = algebra::hypot(dx, dy);

      auto dxsd = dx / d;
      auto dysd = dy / d;

      auto& g1 = h_matrix1;
      auto& g2 = h_matrix2;

      g1(0, 0) = -dxsd;
      g1(1, 0) = dysd;
      g1(0, 1) = -dysd;
      g1(1, 1) = -dxsd;
      g1(0, 2) = 0.;
      g1(1, 2) = -r1();

      g2(0, 0) = dxsd;
      g2(1, 0) = -dysd;
      g2(0, 1) = dysd;
      g2(1, 1) = dxsd;
      g2(0, 2) = 0.;
      g2(1, 2) = -r2();
    }
  };
};

namespace lagrangian {
// free functions for all lagrangian relations

template <typename Data, match::linear_relation HandleRel,
          match::handle<lagrangian_ds> HandleDS>
decltype(auto) compute_h(Data& data, HandleRel& relation, HandleDS& ds)
{
  return relation.h_matrix() * ds.q();
}
template <typename Data, match::item Relation>
decltype(auto) compute_jachq(Data& data, Relation& relation)
{
  if constexpr (has_property_t<Relation, property::time_invariant, Data>{}) {
    return relation.h_matrix();
  }
}
}  // namespace lagrangian

}  // namespace siconos::model
