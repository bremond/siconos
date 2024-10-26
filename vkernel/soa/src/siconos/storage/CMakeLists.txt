add_subdirectory(tests)
add_subdirectory(ground)
add_subdirectory(pattern)
add_subdirectory(some)
add_subdirectory(traits)

add_library(siconos.storage)
target_sources(siconos.storage PUBLIC FILE_SET CXX_MODULES FILES
  siconos_storage.cxx info.cxx add.cxx get.cxx properties.cxx memory.cxx
  make.cxx handle.cxx remove.cxx data_holder.cxx)

target_link_libraries(siconos.storage siconos.storage.ground)
target_link_libraries(siconos.storage siconos.storage.pattern)
target_link_libraries(siconos.storage siconos.storage.some)
target_link_libraries(siconos.storage siconos.storage.traits)

