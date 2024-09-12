#include "py_head.hpp"

void py_add(py::module& mod)
{
  using disks_info_t =
      siconos::storage::get_info_t<siconos::python::disks::idata_t>;

  using disks_properties_t = typename disks_info_t::all_properties_t;

  ground::for_each(py_handles(mod), [&mod](auto pyhandle) {
    using handle_t = typename decltype(+pyhandle[2_c])::type;
    using item_t = typename handle_t::type;
    auto item_name = storage::bind_name<item_t, disks_properties_t>();

    mod.def(
        fmt::format("add_{}", item_name).c_str(),
        [](siconos::python::disks::data_t& data) {
          return siconos::storage::add<item_t>(data());
        },
        py::return_value_policy::reference, R"pbdoc(
        Add
    )pbdoc");
  });
}
