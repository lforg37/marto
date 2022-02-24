# Distributed under the FloPoCo License, see README.md for more information

#[=======================================================================[.rst:
FindSollya
-------

Finds the Sollya library.

Imported Targets
^^^^^^^^^^^^^^^^

This module provides the following imported targets, if found:

``Sollya::Sollya``
  The Sollya library

Result Variables
^^^^^^^^^^^^^^^^

This will define the following variables:

``Sollya_FOUND``
  True if the system has the Sollya library.
``Sollya_VERSION``
  The version of the Sollya library which was found.
``Sollya_INCLUDE_DIRS``
  Include directories needed to use Sollya.
``Sollya_LIBRARIES``
  Libraries needed to link to Sollya.

Cache Variables
^^^^^^^^^^^^^^^

The following cache variables may also be set:

``Sollya_INCLUDE_DIR``
  The directory containing ``sollya.h``.
``Sollya_LIBRARY``
  The path to the Sollya library.

#]=======================================================================]

find_package(PkgConfig)
pkg_check_modules(PC_Sollya QUIET sollya)

find_path(Sollya_INCLUDE_DIR
  NAMES sollya.h
  PATHS ${PC_Sollya_INCLUDE_DIRS}
  DOC "Path of sollya.h, the include file for GNU Sollya library"
)

FIND_LIBRARY(Sollya_LIBRARY
  NAMES libsollya.so libsollya.dylib
  PATHS ${PC_Sollya_LIBRARY_DIRS}
  DOC "Directory of the Sollya library"
)

set(Sollya_VERSION ${PC_Sollya_VERSION})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  Sollya
  FOUND_VAR Sollya_FOUND 
  REQUIRED_VARS
    Sollya_LIBRARY
    Sollya_INCLUDE_DIR
  VERSION_VAR Sollya_VERSION
)

if(Sollya_FOUND)
  set(Sollya_LIBRARIES ${Sollya_LIBRARY})
  set(Sollya_INCLUDE_DIRS ${Sollya_INCLUDE_DIR})
  set(Sollya_DEFINITIONS ${PC_Sollya_FLAGS_OTHER})
endif()

if (Sollya_FOUND AND NOT TARGET Sollya::Sollya)
  find_package(MPFI REQUIRED)
  add_library(Sollya::Sollya UNKNOWN IMPORTED)
  set_target_properties(Sollya::Sollya PROPERTIES 
    IMPORTED_LOCATION "${Sollya_LIBRARY}"
    INTERFACE_COMPILE_OPTIONS "${PC_Sollya_FLAGS_OTHER}"
    INTERFACE_INCLUDE_DIRECTORRIES "${Sollya_INCLUDE_DIR}"
  )
  target_link_libraries(Sollya::Sollya INTERFACE MPFI::MPFI)
endif()

mark_as_advanced(Sollya_INCLUDE_DIR Sollya_LIBRARY)
