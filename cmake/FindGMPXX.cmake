# Distributed under the FloPoCo License, see README.md for more information

#[=======================================================================[.rst:
FindGMPXX
-------

Finds the GMPXX library.

Imported Targets
^^^^^^^^^^^^^^^^

This module provides the following imported targets, if found:

``GMPXX::GMPXX``
  The GMPXX library

Result Variables
^^^^^^^^^^^^^^^^

This will define the following variables:

``GMPXX_FOUND``
  True if the system has the GMPXX library.
``GMPXX_VERSION``
  The version of the GMPXX library which was found.
``GMPXX_INCLUDE_DIRS``
  Include directories needed to use GMPXX.
``GMPXX_LIBRARIES``
  Libraries needed to link to GMPXX.

Cache Variables
^^^^^^^^^^^^^^^

The following cache variables may also be set:

``GMPXX_INCLUDE_DIR``
  The directory containing ``gmpxx.h``.
``GMPXX_LIBRARY``
  The path to the GMPXX library.

#]=======================================================================]

find_package(PkgConfig)
pkg_check_modules(PC_GMPXX QUIET gmpxx)

find_path(GMPXX_INCLUDE_DIR
  NAMES gmpxx.h
  PATHS ${PC_GMPXX_INCLUDE_DIRS}
  DOC "Path of gmp.h, the include file for GNU GMPXX library"
)

FIND_LIBRARY(GMPXX_LIBRARY
  NAMES gmpxx gmpxx.lib
  PATHS ${PC_GMPXX_LIBRARY_DIRS}
  DOC "Directory of the GMPXX library"
)

set(GMPXX_VERSION ${PC_GMPXX_VERSION})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  GMPXX
  FOUND_VAR GMPXX_FOUND 
  REQUIRED_VARS
    GMPXX_LIBRARY
    GMPXX_INCLUDE_DIR
  VERSION_VAR GMPXX_VERSION
)

if(GMPXX_FOUND)
  set(GMPXX_LIBRARIES ${GMPXX_LIBRARY})
  set(GMPXX_INCLUDE_DIRS ${GMPXX_INCLUDE_DIR})
  set(GMPXX_DEFINITIONS ${PC_GMPXX_FLAGS_OTHER})
endif()

if (GMPXX_FOUND AND NOT TARGET GMPXX::GMPXX)
  find_package(GMP REQUIRED)
  add_library(GMPXX::GMPXX UNKNOWN IMPORTED)
  set_target_properties(GMPXX::GMPXX PROPERTIES 
    IMPORTED_LOCATION "${GMPXX_LIBRARY}"
    INTERFACE_COMPILE_OPTIONS "${PC_GMPXX_FLAGS_OTHER}"
    INTERFACE_INCLUDE_DIRECTORRIES "${GMPXX_INCLUDE_DIR}"
  )
  target_link_libraries(GMPXX::GMPXX INTERFACE GMP::GMP)
endif()

mark_as_advanced(GMPXX_INCLUDE_DIR GMPXX_LIBRARY)
