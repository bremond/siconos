#include "py_head.hpp"

void py_methods(py::module& mod)
{
    // methods
  ground::for_each(py_handles(mod), [](auto pyhandle) {
    using handle_t = typename decltype(+pyhandle[2_c])::type;
    //    using item_t = typename handle_t::type;
    ground::fold_left(storage::methods(pyhandle[2_c]),  // all methods
                      std::ref(pyhandle[1_c]),          // initial state
                      []<typename M>(py::class_<handle_t> dc, M m) {
                        return dc.def(pattern::method_name(m),
                                      pattern::method_def(m),
                                      py::return_value_policy::reference);
                      });
  });

}
