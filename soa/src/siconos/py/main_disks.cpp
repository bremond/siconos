#include "py_head.hpp"

void py_data(py::module&);
void py_attributes(py::module&);
void py_attached_storage(py::module&);
void py_methods(py::module&);
void py_add(py::module&);


namespace siconos::python::disks {
data_t make_storage() { return data_t(); };
}  // namespace siconos::python::disks

PYBIND11_MODULE(_nonos, m)
{
  // sub module disks
  auto disks = m.def_submodule("disks");

  py_data(disks);
  py_attributes(disks);
  py_attached_storage(disks);
  py_methods(disks);
  py_add(disks);

  disks.doc() = R"pbdoc(
        Nonos m
        -----------

        .. currentmodule:: nonos

        .. autosummary::
           :toctree: _generate

           add
           subtract
    )pbdoc";

  disks.def("make_storage", &siconos::python::disks::make_storage,
            py::return_value_policy::reference,
            R"pbdoc(
        Create a new data object for 2D disks simulation
    )pbdoc");


  disks.attr("__version__") = "dev";
  m.attr("__version__") = "dev";
}
