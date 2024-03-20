add_executable(t1 t1.cpp)
#add_executable(t2 t2.cpp)
add_executable(bouncing_ball bouncing_ball.cpp)
add_executable(bouncing_balls bouncing_balls.cpp)
add_executable(bouncing_disk bouncing_disk.cpp)
#add_executable(tnonos tnonos.cpp)

set_target_properties(t1 PROPERTIES COMPILE_FLAGS "${flags}")# -fconcepts-diagnostics-depth=10")

#set_target_properties(t2 PROPERTIES COMPILE_FLAGS "${flags}")# -fconcepts-diagnostics-depth=10")
#set_target_properties(tnonos PROPERTIES COMPILE_FLAGS "${flags}")# -fconcepts-diagnostics-depth=10")
set_target_properties(bouncing_ball PROPERTIES COMPILE_FLAGS "${flags}")# -fconcepts-diagnostics-depth=10")
set_target_properties(bouncing_balls PROPERTIES COMPILE_FLAGS "${flags}")# -fconcepts-diagnostics-depth=10")
set_target_properties(bouncing_disk PROPERTIES COMPILE_FLAGS "${flags}")# -fconcepts-diagnostics-depth=10")


target_link_libraries(t1 Eigen3::Eigen Siconos::numerics benchmark pthread fmt)
target_link_libraries(bouncing_ball  Eigen3::Eigen Siconos::numerics benchmark pthread fmt)
target_link_libraries(bouncing_balls Eigen3::Eigen Siconos::numerics benchmark pthread fmt)
target_link_libraries(bouncing_disk Eigen3::Eigen Siconos::numerics benchmark pthread fmt CompactNSearch)



