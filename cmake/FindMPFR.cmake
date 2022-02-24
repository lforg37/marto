# Distributed under the FloPoCo License, see README.md for more information

#[=======================================================================[.rst:
FindMPFR
-------

Finds the MPFR library.

Imported Targets
^^^^^^^^^^^^^^^^

This module provides the following imported targets, if found:

``MPFR::MPFR``
  The MPFR library

Result Variables
^^^^^^^^^^^^^^^^

This will define the following variables:

``MPFR_FOUND``
  True if the system has the MPFR library.
``MPFR_VERSION``
  The version of the MPFR library which was found.
``MPFR_INCLUDE_DIRS``
  Include directories needed to use MPFR.
``MPFR_LIBRARIES``
  Libraries needed to link to MPFR.

Cache Variables
^^^^^^^^^^^^^^^

The following cache variables may also be set:

``MPFR_INCLUDE_DIR``
  The directory containing ``mpfr.h``.
``MPFR_LIBRARY``
  The path to the MPFR library.

#]=======================================================================]

find_package(PkgConfig)
pkg_check_modules(PC_MPFR QUIET mpfr)

find_path(MPFR_INCLUDE_DIR
  NAMES mpfr.h
  PATHS ${PC_MPFR_INCLUDE_DIRS}
  DOC "Path of mpfr.h, the include file for GNU MPFR library"
)

FIND_LIBRARY(MPFR_LIBRARY
  NAMES mpfr mpfr.lib
  PATHS ${PC_MPFR_LIBRARY_DIRS}
  DOC "Directory of the MPFR library"
)

set(MPFR_VERSION ${PC_MPFR_VERSION})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  MPFR
  FOUND_VAR MPFR_FOUND 
  REQUIRED_VARS
    MPFR_LIBRARY
    MPFR_INCLUDE_DIR
  VERSION_VAR MPFR_VERSION
)

if(MPFR_FOUND)
  set(MPFR_LIBRARIES ${MPFR_LIBRARY})
  set(MPFR_INCLUDE_DIRS ${MPFR_INCLUDE_DIR})
  set(MPFR_DEFINITIONS ${PC_MPFR_FLAGS_OTHER})
endif()

if (MPFR_FOUND AND NOT TARGET MPFR::MPFR)
  find_package(GMP REQUIRED)
  add_library(MPFR::MPFR UNKNOWN IMPORTED)
  set_target_properties(MPFR::MPFR PROPERTIES 
    IMPORTED_LOCATION "${MPFR_LIBRARY}"
    INTERFACE_COMPILE_OPTIONS "${PC_MPFR_FLAGS_OTHER}"
    INTERFACE_INCLUDE_DIRECTORRIES "${MPFR_INCLUDE_DIR}"
  )
  target_link_libraries(MPFR::MPFR INTERFACE GMP::GMP)
endif()

mark_as_advanced(MPFR_INCLUDE_DIR MPFR_LIBRARY)
