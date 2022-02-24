# Distributed under the FloPoCo License, see README.md for more information

#[=======================================================================[.rst:
FindGMP
-------

Finds the GMP library.

Imported Targets
^^^^^^^^^^^^^^^^

This module provides the following imported targets, if found:

``GMP::GMP``
  The GMP library

Result Variables
^^^^^^^^^^^^^^^^

This will define the following variables:

``GMP_FOUND``
  True if the system has the GMP library.
``GMP_VERSION``
  The version of the GMP library which was found.
``GMP_INCLUDE_DIRS``
  Include directories needed to use GMP.
``GMP_LIBRARIES``
  Libraries needed to link to GMP.

Cache Variables
^^^^^^^^^^^^^^^

The following cache variables may also be set:

``GMP_INCLUDE_DIR``
  The directory containing ``gmp.h``.
``GMP_LIBRARY``
  The path to the GMP library.

#]=======================================================================]

find_package(PkgConfig)
pkg_check_modules(PC_GMP QUIET gmp)

find_path(GMP_INCLUDE_DIR
  NAMES gmp.h
  PATHS ${PC_GMP_INCLUDE_DIRS}
  DOC "Path of gmp.h, the include file for GNU GMP library"
)

FIND_LIBRARY(GMP_LIBRARY
  NAMES gmp gmp.lib
  PATHS ${PC_GMP_LIBRARY_DIRS}
  DOC "Directory of the GMP library"
)

set(GMP_VERSION ${PC_GMP_VERSION})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  GMP
  FOUND_VAR GMP_FOUND 
  REQUIRED_VARS
    GMP_LIBRARY
    GMP_INCLUDE_DIR
  VERSION_VAR GMP_VERSION
)

if(GMP_FOUND)
  set(GMP_LIBRARIES ${GMP_LIBRARY})
  set(GMP_INCLUDE_DIRS ${GMP_INCLUDE_DIR})
  set(GMP_DEFINITIONS ${PC_GMP_FLAGS_OTHER})
endif()

if (GMP_FOUND AND NOT TARGET GMP::GMP)
  add_library(GMP::GMP UNKNOWN IMPORTED)
  set_target_properties(GMP::GMP PROPERTIES 
    IMPORTED_LOCATION "${GMP_LIBRARY}"
    INTERFACE_COMPILE_OPTIONS "${PC_GMP_FLAGS_OTHER}"
    INTERFACE_INCLUDE_DIRECTORRIES "${GMP_INCLUDE_DIR}"
  )
endif()

mark_as_advanced(GMP_INCLUDE_DIR GMP_LIBRARY)
