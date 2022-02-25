cmake_minimum_required(VERSION 3.10)
project(Marto CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Add module to find custom dependencies
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

include(FetchContent)

set(LIBNAME marto)
set(INCLUDE_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/include)

if (NOT DEFINED HintBackendVersion)
	set(HintBackendVersion "master")
endif()
set (HintBackendVersion ${HintBackendVersion} CACHE STRING "Hint commit version used")

FetchContent_Declare(
	Hint
	GIT_REPOSITORY https://github.com/lforg37/hint.git
	GIT_TAG ${HintBackendVersion}
	)

FetchContent_MakeAvailable(Hint)

add_library(${LIBNAME} INTERFACE)
target_include_directories(${LIBNAME} INTERFACE 
	$<BUILD_INTERFACE:${INCLUDE_ROOT}>
	$<INSTALL_INTERFACE:include/marto>
	)
target_link_libraries(${LIBNAME} INTERFACE hint)

add_subdirectory(libnumform)

OPTION(BUILD_ARCHGENLIB "Build archgenlib" ON)
if(BUILD_ARCHGENLIB)
add_subdirectory(archgenlib)
endif()


OPTION(DEBUG_FPEXPR_SUM "Dump all FPExpr sum operator signals")
if(DEBUG_FPEXPR_SUM)
	target_compile_definitions(${LIBNAME} INTERFACE FPEXPR_SUM_DEBUG)
endif()

OPTION(DEBUG_FPEXPR_ROUND "Dump all FPExpr Round operator signals")
if(DEBUG_FPEXPR_ROUND)
	target_compile_definitions(${LIBNAME} INTERFACE FPEXPR_ROUND_DEBUG)
endif()

OPTION(DEBUG_IEEE_ADDER "Dump all ieee_adder signals")
if(DEBUG_IEEE_ADDER)
	target_compile_definitions(${LIBNAME} INTERFACE IEEE_ADDER_DEBUG)
endif()

OPTION(DEBUG_IEEE_ADDER_SORTEXP "Dump all ieee_adder_sortexp signals")
if(DEBUG_IEEE_ADDER_SORTEXP)
	target_compile_definitions(${LIBNAME} INTERFACE IEEE_ADDER_SORTEXP_DEBUG)
endif()

OPTION(DEBUG_POSIT_ADDER "Dump all posit_adder signals")
if(DEBUG_POSIT_ADDER)
	target_compile_definitions(${LIBNAME} INTERFACE POSIT_ADDER_DEBUG)
endif()

OPTION(DEBUG_POSIT_DECODER "Dump all posit_decoder signals")
if(DEBUG_POSIT_DECODER)
	target_compile_definitions(${LIBNAME} INTERFACE POSIT_DECODER_DEBUG)
endif()

OPTION(DEBUG_POSIT_ENCODER "Dump all posit_encoder signals")
if(DEBUG_POSIT_ENCODER)
	target_compile_definitions(${LIBNAME} INTERFACE POSIT_ENCODER_DEBUG)
endif()

OPTION(DEBUG_POSIT_MUL "Dump all posit_mul signals")
if(DEBUG_POSIT_MUL)
	target_compile_definitions(${LIBNAME} INTERFACE POSIT_MUL_DEBUG)
endif()

OPTION(DEBUG_POSIT_QUIRE_ADD "Dump all posit  add_sub_quire signals")
if(DEBUG_POSIT_QUIRE_ADD)
	target_compile_definitions(${LIBNAME} INTERFACE POSIT_QUIRE_ADD_DEBUG)
endif()

OPTION(DEBUG_POSIT_QUIRE_PIF "Dump all posit  quire_to_posit signals")
if(DEBUG_POSIT_QUIRE_PIF)
	target_compile_definitions(${LIBNAME} INTERFACE POSIT_QUIRETOPIF_DEBUG)
endif()

OPTION(DEBUG_POSIT_ROUNDER "Dump all posit  posit_round_in_place signals")
if(DEBUG_POSIT_ROUNDER)
	target_compile_definitions(${LIBNAME} INTERFACE POSIT_ROUNDER_DEBUG)
endif()

OPTION(DEBUG_POSIT_VALPROD "Dump all value_prod_conversion signals")
if(DEBUG_POSIT_VALPROD)
	target_compile_definitions(${LIBNAME} INTERFACE POSIT_VALPROD_DEBUG)
endif()

OPTION(DEBUG_IEEE_MUL_exact "Dump all fp_exact_prod signals")
if(DEBUG_IEEE_MUL_exact)
	target_compile_definitions(${LIBNAME} INTERFACE IEEE_MUL_DEBUG)
endif()

OPTION(DEBUG_K1_ACC "Dump all k1 accumulation signals")
if(DEBUG_K1_ACC)
	target_compile_definitions(${LIBNAME} INTERFACE K1_ACC_DEBUG)
endif()

OPTION(DEBUG_K3_ACC "Dump all k3 accumulation signals")
if(DEBUG_K3_ACC)
	target_compile_definitions(${LIBNAME} INTERFACE K3_ACC_DEBUG)
endif()

OPTION(DEBUG_KULISH_ROUND "Dump all kulish rounding signals")
if(DEBUG_KULISH_ROUND)
	target_compile_definitions(${LIBNAME} INTERFACE KULISCH_ROUND_DEBUG)
endif()

install(DIRECTORY ${INCLUDE_ROOT}/ DESTINATION include/marto FILES_MATCHING PATTERN "*.hpp" PATTERN "*.ipp")
install(
	TARGETS ${LIBNAME}
	EXPORT ${PROJECT_NAME}Targets
	PUBLIC_HEADER DESTINATION "include/marto"
	)
install(FILES ${PROJECT_NAME}Config.cmake DESTINATION lib/cmake/${PROJECT_NAME})
install(EXPORT ${PROJECT_NAME}Targets FILE ${PROJECT_NAME}Targets.cmake DESTINATION lib/cmake/${PROJECT_NAME})
export(PACKAGE ${PROJECT_NAME})

OPTION(BUILD_MARTO_UNIT_TESTS "Build the operators unit tests")
if (BUILD_MARTO_UNIT_TESTS)
    add_subdirectory(tests)
endif()
