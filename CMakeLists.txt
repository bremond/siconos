cmake_minimum_required(VERSION 3.2)
project(mybenchmarks_library)
include(CheckCXXCompilerFlag)
#add_compile_options(-fPIC)
#add_compile_options(-fPIE)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
include_directories(${CMAKE_SOURCE_DIR}/../ctti/include)
include_directories(${CMAKE_SOURCE_DIR}/numeric_bindings)
include_directories(
  ${SICONOS_INSTALL_DIR}/include
  ${SICONOS_INSTALL_DIR}/include/siconos
  ${SICONOS_INSTALL_DIR}/include/siconos/numerics
  ${SICONOS_INSTALL_DIR}/include/siconos/kernel)
link_directories(${CMAKE_BUILD_DIR})
add_library(xsiconos SHARED XSiconosVector.cpp)
add_executable(tbench tbench.cpp)
#add_executable(tgraph tgraph.cpp)

set(flags "-Wall --std=c++20")

set(fconcepts_diagnostics)
check_cxx_compiler_flag("-fconcepts-diagnostics-depth=0" fconcepts_diagnostics)

if(fconcepts_diagnostics)
  set(flags "${flags} -fconcepts-diagnostics-depth=0")
endif()
#add_executable(bouncing_ball bouncing_ball.cpp)
set_target_properties(xsiconos PROPERTIES COMPILE_FLAGS "${flags}")
set_target_properties(tbench PROPERTIES COMPILE_FLAGS "${flags}")
#set_target_properties(tgraph PROPERTIES COMPILE_FLAGS "${flags}")

#set_target_properties(bouncing_ball PROPERTIES COMPILE_FLAGS "-Wall --std=c++20")
target_link_libraries(xsiconos siconos_kernel fmt)
target_link_libraries(tbench xsiconos benchmark pthread fmt)
#target_link_libraries(tgraph xsiconos benchmark pthread fmt)
target_link_libraries(tnonos xsiconos benchmark pthread fmt)
#target_link_libraries(bouncing_ball xsiconos benchmark pthread fmt)


