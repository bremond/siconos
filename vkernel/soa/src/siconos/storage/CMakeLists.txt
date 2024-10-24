add_subdirectory(tests)

add_library(siconos.storage)
target_sources(siconos.storage PUBLIC FILE_SET CXX_MODULES FILES
  siconos_storage.cxx info.cxx get.cxx properties.cxx)
