#include "py_head.hpp"

static_assert(storage::pattern::match::diagonal_matrix<
              decltype(storage::get<
                       storage::pattern::attr_t<config::disk, "mass_matrix">>(
                  idata_t{}, 0, storage::add<config::disk>(idata_t{})))>);


void py_data(py::module mod&)
{
  auto data_class =
      py::class_<siconos::python::disks::data_t>(mod, "data_t");
}
