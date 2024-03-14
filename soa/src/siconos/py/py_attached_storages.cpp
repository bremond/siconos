#include "py_head.hpp"

void py_attached_storages(py::module& mod)
{
  // attached storages
  ground::for_each(py_handles(mod), [](auto pyhandle) {
    using handle_t = typename decltype(+pyhandle[2_c])::type;
    using item_t = typename handle_t::type;

    ground::fold_left(
        decltype(storage::attached_storages(
            handle_t{}, handle_t{}.data())){},  // all attached storages
        std::ref(pyhandle[1_c]),                // initial state
        []<match::attached_storage<item_t> S>(py::class_<handle_t> dc, S s) {
          constexpr auto astor_name = storage::attached_storage_name(s);
          using target_type = std::decay_t<decltype(out_formatter(
              handle_t{}, storage::prop<astor_name.str>(handle_t{})))>;

          return dc
              .def(
                  fmt::format("{}", astor_name.str.value).c_str(),
                  [&astor_name](handle_t& h) -> target_type {
                    return out_formatter(h, storage::prop<astor_name.str>(h));
                  },
                  py::return_value_policy::reference)
              .def(fmt::format("set_{}", astor_name.str.value).c_str(),
                   [&astor_name](handle_t& h, target_type val) {
                     in_formatter(h, storage::prop<astor_name.str>(h)) =
                         in_formatter(h, val);
                   });
        });
  });
}
