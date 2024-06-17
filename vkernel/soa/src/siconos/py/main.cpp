#include "py_head.hpp"

namespace siconos::python::disks {
data_t make_storage() { return data_t(); };
}  // namespace siconos::python::disks

using namespace boost::hana::literals;
PYBIND11_MODULE(_nonos, m)
{
  // sub module disks
  auto disks = m.def_submodule("disks");

  auto data_class =
      py::class_<siconos::python::disks::data_t>(disks, "data_t");

  using disks_info_t = std::decay_t<decltype(ground::get<storage::info>(
      siconos::python::disks::idata_t{}))>;

  using disks_properties_t = typename disks_info_t::all_properties_t;

  using disks_items_t = decltype(ground::transform(
      typename disks_info_t::all_items_t{}, []<match::item I>(I) {
        if constexpr (match::wrap<I>) {
          return typename I::type{};
        }
        else {
          return I{};
        }
      }));

  // ground::type_trace<disks_items_t>();
  auto named_disks_items = ground::tuple_unique(ground::filter(
      disks_items_t{}, ground::is_a_model<[]<typename T>() {
        return storage::has_property_from<T, storage::property::bind,
                                          disks_properties_t>();
      }>));

  // ground::type_trace<std::decay_t<decltype(named_disks_items)>>();

  auto disks_handles = ground::transform(
      // only named items
      named_disks_items, []<match::item I>(I item) {
        // get handle
        return storage::add<I>(siconos::python::disks::idata_t{});
      });

  // add corresponding py::class_
  auto pyhandles = ground::transform(disks_handles, [&disks]<typename H>(
                                                        H handle) {
    using item_t = typename H::type;
    using base_index_t = typename H::base_index_t;
    auto base_index = py::class_<base_index_t>(
        disks, fmt::format("index_{}",
                           storage::bind_name<item_t, disks_properties_t>())
                   .c_str());
    return ground::make_tuple(
        base_index,
        py::class_<H>(disks, storage::bind_name<item_t, disks_properties_t>(),
                      base_index),
        ground::type_c<H>);
  });

  // attached storage
  ground::for_each(pyhandles, [](auto pyhandle) {
    using handle_t = typename decltype(+pyhandle[2_c])::type;
    using item_t = typename handle_t::type;

    ground::fold_left(
        decltype(storage::attached_storages(
            item_t{}, handle_t{}.data())){},  // all attached storages
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

  // attributes
  ground::for_each(pyhandles, [](auto pyhandle) {
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

  // methods
  ground::for_each(pyhandles, [](auto pyhandle) {
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

  ground::for_each(pyhandles, [&disks](auto pyhandle) {
    using handle_t = typename decltype(+pyhandle[2_c])::type;
    using item_t = typename handle_t::type;
    auto item_name = storage::bind_name<item_t, disks_properties_t>();

    disks.def(
        fmt::format("add_{}", item_name).c_str(),
        [](siconos::python::disks::data_t& data) {
          return siconos::storage::add<item_t>(data());
        },
        py::return_value_policy::reference, R"pbdoc(
        Add
    )pbdoc");
  });

  disks.attr("__version__") = "dev";
  m.attr("__version__") = "dev";
}
