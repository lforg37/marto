# Distributed under the FloPoCo License, see README.md for more information

#[=======================================================================[.rst:
FindMPFI
-------

Finds the MPFI library.

Imported Targets
^^^^^^^^^^^^^^^^

This module provides the following imported targets, if found:

``MPFI::MPFI``
  The MPFI library

Result Variables
^^^^^^^^^^^^^^^^

This will define the following variables:

``MPFI_FOUND``
  True if the system has the MPFI library.
``MPFI_VERSION``
  The version of the MPFI library which was found.
``MPFI_INCLUDE_DIRS``
  Include directories needed to use MPFI.
``MPFI_LIBRARIES``
  Libraries needed to link to MPFI.

Cache Variables
^^^^^^^^^^^^^^^

The following cache variables may also be set:

``MPFI_INCLUDE_DIR``
  The directory containing ``mpfi.h``.
``MPFI_LIBRARY``
  The path to the MPFI library.

#]=======================================================================]

find_package(PkgConfig)
pkg_check_modules(PC_MPFI QUIET mpfi)

find_path(MPFI_INCLUDE_DIR
  NAMES mpfi.h
  PATHS ${PC_MPFI_INCLUDE_DIRS}
  DOC "Path of mpfi.h, the include file for GNU MPFI library"
)

FIND_LIBRARY(MPFI_LIBRARY
  NAMES mpfi mpfi.lib
  PATHS ${PC_MPFI_LIBRARY_DIRS}
  DOC "Directory of the MPFI library"
)

set(MPFI_VERSION ${PC_MPFI_VERSION})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  MPFI
  FOUND_VAR MPFI_FOUND 
  REQUIRED_VARS
    MPFI_LIBRARY
    MPFI_INCLUDE_DIR
  VERSION_VAR MPFI_VERSION
)

if(MPFI_FOUND)
  set(MPFI_LIBRARIES ${MPFI_LIBRARY})
  set(MPFI_INCLUDE_DIRS ${MPFI_INCLUDE_DIR})
  set(MPFI_DEFINITIONS ${PC_MPFI_FLAGS_OTHER})
endif()

if (MPFI_FOUND AND NOT TARGET MPFI::MPFI)
  find_package(GMP REQUIRED)
  find_package(MPFR REQUIRED)
  add_library(MPFI::MPFI UNKNOWN IMPORTED)
  set_target_properties(MPFI::MPFI PROPERTIES 
    IMPORTED_LOCATION "${MPFI_LIBRARY}"
    INTERFACE_COMPILE_OPTIONS "${PC_MPFI_FLAGS_OTHER}"
    INTERFACE_INCLUDE_DIRECTORRIES "${MPFI_INCLUDE_DIR}"
  )
  target_link_libraries(MPFI::MPFI INTERFACE GMP::GMP MPFI::MPFI)
endif()

mark_as_advanced(MPFI_INCLUDE_DIR MPFI_LIBRARY)
