# Find MUMPS includes and libraries.
# The following variables are set if MUMPS is found.  If MUMPS is not
# found, MUMPS_FOUND is set to false.
#  MUMPS_FOUND        - True when the MUMPS include directory is found.
#  MUMPS_INCLUDE_DIRS - the path to where the Siconos MUMPS include files are.
#  MUMPS_LIBRARY_DIRS - The path to where the Siconos library files are.
#  MUMPS_LIBRARIES    - The libraries to link against Siconos MUMPS

# One may want to use a specific MUMPS Library by setting
# MUMPS_LIBRARY_DIRECTORY before FIND_PACKAGE(MUMPS)
INCLUDE(FindPackageHandleStandardArgs)

# On debian and Suse, both mpi and mpi-free version of MUMPS can be installed in parallel
# The latter have a "_seq" suffix
IF(IDONTWANTMPI)
  SET(__MUMPS_NAMES dmumps_seq dmumps)
ELSE(IDONTWANTMPI)
# on debian systems one may have mumps+[pt]scotch packages also
  SET(__MUMPS_NAMES dmumps_ptscotch dmumps_scotch dmumps)
ENDIF(IDONTWANTMPI)

IF(MUMPS_LIBRARY_DIRECTORY)
  MESSAGE(STATUS "Searching for ${__MUMPS_NAMES} in ${MUMPS_LIBRARY_DIRECTORY}")
  FIND_LIBRARY(MUMPS_LIBRARY NAMES ${__MUMPS_NAMES} PATHS "${MUMPS_LIBRARY_DIRECTORY}" NO_DEFAULT_PATH)
  IF(NOT MUMPS_LIBRARY)
    SET(_MUMPS_DEFAULT_SEARCH_DIR TRUE)
    MESSAGE("Could not find MUMPS in directory ${MUMPS_LIBRARY_DIRECTORY}")
  ENDIF()
ENDIF()

IF(NOT MUMPS_LIBRARY)
  MESSAGE(STATUS "Searching for ${__MUMPS_NAMES} in default directories")
  FIND_LIBRARY(MUMPS_LIBRARY NAMES ${__MUMPS_NAMES})
ENDIF()

IF(MUMPS_LIBRARY)
  MESSAGE(STATUS "Found ${MUMPS_LIBRARY}")
ENDIF()

# Try to be smart and detect whether the MUMPS lib file has a "_seq" suffix. If yes, we add it to all the other library
IF(MUMPS_LIBRARY MATCHES "dmumps_seq")
  SET(__MUMPS_COMMON_NAMES mumps_common_seq mumps_common)
ELSE()
  SET(__MUMPS_COMMON_NAMES mumps_common_ptscotch mumps_common_scotch mumps_common mumps_common_seq)
ENDIF()

IF(MUMPS_LIBRARY_DIRECTORY AND NOT _MUMPS_DEFAULT_SEARCH_DIR)
  MESSAGE(STATUS "Searching for ${__MUMPS_COMMON_NAMES} in ${MUMPS_LIBRARY_DIRECTORY}")
  FIND_LIBRARY(MUMPS_COMMON_LIBRARY NAMES ${__MUMPS_COMMON_NAMES} PATHS "${MUMPS_LIBRARY_DIRECTORY}" NO_DEFAULT_PATH)
ENDIF()

IF(NOT MUMPS_COMMON_LIBRARY)
  MESSAGE(STATUS "Searching for ${__MUMPS_COMMON_NAMES} in default directories")
  FIND_LIBRARY(MUMPS_COMMON_LIBRARY NAMES ${__MUMPS_COMMON_NAMES})
ENDIF()

IF(MUMPS_COMMON_LIBRARY)
  MESSAGE(STATUS "Found ${MUMPS_COMMON_LIBRARY}")
ENDIF()

IF(MUMPS_LIBRARY_DIRECTORY)
  FIND_LIBRARY(PTSCOTCH_LIBRARY ptscotch PATHS "${MUMPS_LIBRARY_DIRECTORY}" NO_DEFAULT_PATH)
  IF(NOT PTSCOTCH_LIBRARY)
    FIND_LIBRARY(PTSCOTCH_LIBRARY ptscotch)
  ENDIF()
ELSE(MUMPS_LIBRARY_DIRECTORY)
  FIND_LIBRARY(PTSCOTCH_LIBRARY ptscotch)
ENDIF(MUMPS_LIBRARY_DIRECTORY)

IF(MUMPS_LIBRARY_DIRECTORY)
  FIND_LIBRARY(SCOTCH_LIBRARY scotch PATHS "${MUMPS_LIBRARY_DIRECTORY}" NO_DEFAULT_PATH)
  IF(NOT SCOTCH_LIBRARY)
    FIND_LIBRARY(SCOTCH_LIBRARY scotch)
  ENDIF()
ELSE(MUMPS_LIBRARY_DIRECTORY)
  FIND_LIBRARY(SCOTCH_LIBRARY scotch)
ENDIF(MUMPS_LIBRARY_DIRECTORY)


IF(MUMPS_LIBRARY_DIRECTORY)
  FIND_LIBRARY(METIS_LIBRARY metis PATHS "${MUMPS_LIBRARY_DIRECTORY}" NO_DEFAULT_PATH)
  IF(NOT METIS_LIBRARY)
    FIND_LIBRARY(METIS_LIBRARY metis)
  ENDIF()
ELSE(MUMPS_LIBRARY_DIRECTORY)
  FIND_LIBRARY(METIS_LIBRARY metis)
ENDIF(MUMPS_LIBRARY_DIRECTORY)

IF(MUMPS_LIBRARY_DIRECTORY)
  FIND_LIBRARY(PORD_LIBRARY pord${__SUFFIX} PATHS "${MUMPS_LIBRARY_DIRECTORY}" NO_DEFAULT_PATH)
  IF(NOT PORD_LIBRARY)
    FIND_LIBRARY(PORD_LIBRARY pord${__SUFFIX})
  ENDIF()
ELSE(MUMPS_LIBRARY_DIRECTORY)
  FIND_LIBRARY(PORD_LIBRARY pord${__SUFFIX})
ENDIF(MUMPS_LIBRARY_DIRECTORY)

IF((NOT MUMPS_INCLUDE_DIR) AND MUMPS_LIBRARY)
  GET_FILENAME_COMPONENT(MUMPS_LIBRARY_DIR ${MUMPS_LIBRARY} PATH)
  GET_FILENAME_COMPONENT(MUMPS_LIBRARY_DIR_DIR ${MUMPS_LIBRARY_DIR} PATH)

  # Suse uses /usr/include/mumps
  # Debian has /usr/include/mumps_seq for the MPI-free version
  find_path(MUMPS_INCLUDE_DIR dmumps_c.h
    HINTS
    ${MUMPS_LIBRARY_DIR_DIR}/include
    $ENV{MPI_INCLUDE}
    PATH_SUFFIXES MUMPS mumps mumps${__SUFFIX}
    )
ENDIF()


FIND_PACKAGE_HANDLE_STANDARD_ARGS(MUMPS
  REQUIRED_VARS MUMPS_LIBRARY MUMPS_INCLUDE_DIR)

IF(MUMPS_LIBRARY)
  GET_FILENAME_COMPONENT(MUMPS_LIBRARY_DIRS ${MUMPS_LIBRARY} PATH)
  SET(MUMPS_LIBRARIES ${MUMPS_LIBRARY} ${MUMPS_COMMON_LIBRARY})

  IF(PTSCOTCH_LIBRARY)
    SET(MUMPS_LIBRARIES ${MUMPS_LIBRARIES} ${PTSCOTCH_LIBRARY})
  ENDIF(PTSCOTCH_LIBRARY)

  IF(SCOTCH_LIBRARY)
    SET(MUMPS_LIBRARIES ${MUMPS_LIBRARIES} ${SCOTCH_LIBRARY})
  ENDIF(SCOTCH_LIBRARY)

  IF(METIS_LIBRARY)
    SET(MUMPS_LIBRARIES ${MUMPS_LIBRARIES} ${METIS_LIBRARY})
  ENDIF(METIS_LIBRARY)

  IF(PORD_LIBRARY)
    SET(MUMPS_LIBRARIES ${MUMPS_LIBRARIES} ${PORD_LIBRARY})
  ENDIF()

  SET(MUMPS_INCLUDE_DIRS ${MUMPS_INCLUDE_DIR})

ELSE(MUMPS_LIBRARY)
  IF(MUMPS_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR
      "Required MUMPS library not found. Please specify library location in MUMPS_LIBRARY_DIRECTORY")
  ENDIF(MUMPS_FIND_REQUIRED)
ENDIF(MUMPS_LIBRARY)
