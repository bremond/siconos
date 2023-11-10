
python_add_library(_nonos MODULE main.cpp WITH_SOABI)

set_target_properties(_nonos PROPERTIES COMPILE_FLAGS "${flags}")
target_link_libraries(_nonos PRIVATE Eigen3::Eigen Siconos::numerics benchmark pthread fmt pybind11::module)
target_compile_definitions(_nonos PRIVATE VERSION_INFO=${PROJECT_VERSION})

install(TARGETS _nonos DESTINATION nonos)

