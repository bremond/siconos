cmake_minimum_required(VERSION 3.2)
project(
  ${SKBUILD_PROJECT_NAME}
  VERSION ${SKBUILD_PROJECT_VERSION}
  LANGUAGES CXX C)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(CheckCXXCompilerFlag)
set(Siconos_DIR ${SICONOS_INSTALL_DIR})

find_package(Siconos REQUIRED)
find_package(Boost REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)

include_directories(${CMAKE_SOURCE_DIR}/src)
link_directories(${CMAKE_BUILD_DIR})

set(flags "-Wall --std=c++20")
set(fconcepts_diagnostics)
check_cxx_compiler_flag("-fconcepts-diagnostics-depth=0" fconcepts_diagnostics)
check_cxx_compiler_flag("-ftemplate-backtrace-limit=0" ftemplate_backtrace_limit)
if(fconcepts_diagnostics)
  set(flags "${flags} -fconcepts-diagnostics-depth=1")
endif()
if(ftemplate_backtrace_limit)
  set(flags "${flags} -ftemplate-backtrace-limit=0")
endif()

# python_add_library(_nothing MODULE ./src/siconos/pybind/nothing.cpp WITH_SOABI)
# set_target_properties(_nothing PROPERTIES COMPILE_FLAGS "${flags}")
# target_link_libraries(_nothing PRIVATE Eigen3::Eigen Siconos::numerics benchmark pthread fmt pybind11::headers)
# target_compile_definitions(_nothing PRIVATE VERSION_INFO=${PROJECT_VERSION})

# install(TARGETS _nothing DESTINATION nothing)

python_add_library(_nonos MODULE ./src/siconos/pybind/main.cpp WITH_SOABI)
set_target_properties(_nonos PROPERTIES COMPILE_FLAGS "${flags}")
target_link_libraries(_nonos PRIVATE Eigen3::Eigen Siconos::numerics benchmark pthread fmt pybind11::headers)
target_compile_definitions(_nonos PRIVATE VERSION_INFO=${PROJECT_VERSION})

install(TARGETS _nonos DESTINATION nonos)

add_subdirectory(src)

