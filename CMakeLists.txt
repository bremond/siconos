cmake_minimum_required(VERSION 3.2)
project(mybenchmarks_library)
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
add_executable(tgraph tgraph.cpp)
add_executable(tnonos tnonos.cpp)
#add_executable(bouncing_ball bouncing_ball.cpp)
set_target_properties(xsiconos PROPERTIES COMPILE_FLAGS "-Wall --std=c++20")
set_target_properties(tbench PROPERTIES COMPILE_FLAGS "-Wall --std=c++20")
set_target_properties(tgraph PROPERTIES COMPILE_FLAGS "-Wall --std=c++20")
set_target_properties(tnonos PROPERTIES COMPILE_FLAGS "-Wall --std=c++20 -fconcepts-diagnostics-depth=10")
#set_target_properties(bouncing_ball PROPERTIES COMPILE_FLAGS "-Wall --std=c++20")
target_link_libraries(xsiconos siconos_kernel fmt)
target_link_libraries(tbench xsiconos benchmark pthread fmt)
target_link_libraries(tgraph xsiconos benchmark pthread fmt)
target_link_libraries(tnonos xsiconos benchmark pthread fmt)
#target_link_libraries(bouncing_ball xsiconos benchmark pthread fmt)


