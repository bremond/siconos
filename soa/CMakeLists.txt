cmake_minimum_required(VERSION 3.2)
project(siconos)

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
  set(flags "${flags} -fconcepts-diagnostics-depth=0")
endif()
if(ftemplate_backtrace_limit)
  set(flags "${flags} -ftemplate-backtrace-limit=0")
endif()

add_subdirectory(src)

