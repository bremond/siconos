#include "py_head.hpp"

void py_attributes(py::module& mod)
{
    // attributes
  ground::for_each(py_handles(mod), [](auto pyhandle) {
    using handle_t = typename decltype(+pyhandle[2_c])::type;
    //    using item_t = typename handle_t::type;
    ground::fold_left(
        pattern::attributes(typename handle_t::type{}),  // all attributes
        std::ref(pyhandle[1_c]),                         // initial state
        []<match::attribute A>(py::class_<handle_t> dc, A a) {
          using attr_value_t = std::decay_t<decltype(out_formatter(
              handle_t{},
              storage::get<A>(handle_t{}.data(), 0, handle_t{})))>;
          return dc
              .def(
                  fmt::format("{}", pattern::attribute_name(a)).c_str(),
                  [](handle_t& h) -> attr_value_t {
                    return out_formatter(h, storage::get<A>(h.data(), h));
                  },
                  py::return_value_policy::reference)
              .def(fmt::format("set_{}", pattern::attribute_name(a)).c_str(),
                   [](handle_t& h, attr_value_t v) {
                     in_formatter(h, storage::get<A>(h.data(), h)) =
                         in_formatter(h, v);
                   });
        });
  });

}
