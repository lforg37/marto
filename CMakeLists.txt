cmake_minimum_required(VERSION 3.10)
project(Marto CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(LIBNAME marto)
set(INCLUDE_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/include)

find_package(libhint CONFIG REQUIRED)

if(NOT DEFINED AP_INT_INCLUDE_DIR)
	MESSAGE(FATAL_ERROR "Please provide the path to ap_int include directory in AP_INT_INCLUDE_DIR")
endif()

set(AP_INT_INCLUDE_DIR ${AP_INT_INCLUDE_DIR} CACHE PATH "Path to ap_int include directory" FORCE)

add_library(VIVADO_AP_INT INTERFACE IMPORTED)
target_include_directories(VIVADO_AP_INT INTERFACE ${AP_INT_INCLUDE_DIR})

add_library(${LIBNAME} INTERFACE)
target_include_directories(${LIBNAME} INTERFACE ${INCLUDE_ROOT})
target_link_libraries(${LIBNAME} INTERFACE hint)


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

OPTION(BUILD_UNIT_TEST "Build the operators unit tests")
if(BUILD_UNIT_TEST)
	find_package(Boost 1.55 REQUIRED COMPONENTS unit_test_framework)
	find_package(OpenMP)

	enable_testing()

	set(TESTDATAROOT ${CMAKE_CURRENT_SOURCE_DIR}/testdata)

	add_executable(testConvert_exe tests/testConvert.cpp  tests/testConversionPosit.cpp)
	find_library(SOFTPOSIT_LIB softposit)
	find_path(SOFTPOSIT_H softposit.h)
	target_link_libraries(testConvert_exe ${LIBNAME} ${Boost_LIBRARIES})
	add_dependencies(testConvert_exe ${LIBNAME})
	if(SOFTPOSIT_LIB AND SOFTPOSIT_H)
		MESSAGE(STATUS "softposit found : ${SOFTPOSIT_LIB}, ${SOFTPOSIT_H}")
		target_compile_definitions(testConvert_exe PRIVATE SOFTPOSIT)
		target_include_directories(testConvert_exe PRIVATE ${SOFTPOSIT_H})
		target_link_libraries(testConvert_exe ${SOFTPOSIT_LIB} VIVADO_AP_INT)
		add_test(testConvert testConvert_exe)

		add_executable(testAritOps_exe tests/testArithmeticOp.cpp)
		target_include_directories(testAritOps_exe PRIVATE ${SOFTPOSIT_H})
		target_link_libraries(testAritOps_exe PRIVATE ${LIBNAME} ${Boost_LIBRARIES} ${SOFTPOSIT_LIB} VIVADO_AP_INT -fopenmp)
		target_compile_definitions(testAritOps_exe PRIVATE TESTDATA_ROOT="${TESTDATAROOT}" VIVADO_BACKEND)
		add_test(testAritOps testAritOps_exe)

		if(OpenMP_CXX_FOUND)
			target_link_libraries(testAritOps_exe PRIVATE OpenMP::OpenMP_CXX)
		endif()
	else()
		Message(WARNING "Tests on posit value cannot be run")
	endif()

	add_executable(testExpressions_exe tests/testExpression.cpp tests/testExprSum.cpp tests/testExprProduct.cpp tests/testExprRounding.cpp)
	target_link_libraries(testExpressions_exe ${LIBNAME} ${Boost_LIBRARIES} VIVADO_AP_INT)
	target_compile_definitions(testExpressions_exe PRIVATE VIVADO_BACKEND)

	find_path(MPFR mpfr.h)

	if(MPFR)
		MESSAGE(STATUS "MPFR found : ${MPFR}")
		add_executable(testKulisch_exe tests/testKulisch.cpp)
		target_link_libraries(testKulisch_exe PRIVATE ${LIBNAME} ${Boost_LIBRARIES} -lmpfr)
		target_compile_definitions(testKulisch_exe PRIVATE MPFR VIVADO_BACKEND)
		target_link_libraries(testExpressions_exe -lmpfr)
		target_compile_definitions(testExpressions_exe PRIVATE MPFR)
		add_test(testKulisch testKulisch_exe)
	else()
		Message(WARNING "Tests on Kulisch value cannot be run")
	endif()

	find_path(SOFTFLOAT_INCLUDE_DIR softfloat.h PATH_SUFFIXES softfloat)
	find_library(SOFTFLOAT libsoftfloat-3.a)

	if (NOT (SOFTFLOAT_INCLUDE_DIR-NOTFOUND OR SOFTFLOAT-NOTFOUND))
		set(SOFTFLOATTARGET-OK ON)
		message(STATUS "Softfloat found")
		message(STATUS "Softfloat lib : ${SOFTFLOAT}")
	else()
		set(SOFTFLOATTARGET-NOTFOUND)
		message(warning "Softfloat not found")
	endif()

	add_library(libnumform INTERFACE)
	target_include_directories(libnumform INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/libnumform)

	add_executable(testLibNumForm_exe tests/testLibNumForm.cpp)
	target_link_libraries(testLibNumForm_exe PRIVATE libnumform ${Boost_LIBRARIES})
	add_test(testLibNumForm testLibNumForm_exe)

	add_executable(ieee_adder_exe tests/ieee_adder.cpp)
	target_link_libraries(ieee_adder_exe PRIVATE VIVADO_AP_INT ${LIBNAME})
	target_compile_definitions(ieee_adder_exe PRIVATE -DIEEE_WE=8 -DIEEE_WF=23)

	add_executable(testIEEElibnum_exe tests/testIEEEMultiplier.cpp tests/testIEEEAdder.cpp)
	target_compile_definitions(testIEEElibnum_exe PRIVATE VIVADO_BACKEND)
	target_link_libraries(testIEEElibnum_exe PRIVATE ${LIBNAME} libnumform ${Boost_LIBRARIES} VIVADO_AP_INT)
	add_test(testIEEElibnum testIEEElibnum_exe)

	if(OpenMP_CXX_FOUND)
		target_link_libraries(testIEEElibnum_exe PRIVATE OpenMP::OpenMP_CXX)
	endif()

	if(SOFTFLOATTARGET-OK)
		target_compile_definitions(testIEEElibnum_exe PRIVATE SOFTFLOAT)
		target_link_libraries(testIEEElibnum_exe PRIVATE ${SOFTFLOAT})
		target_include_directories(testIEEElibnum_exe PRIVATE ${SOFTFLOAT_INCLUDE_DIR})
	else()
		message(warning "Softfloat not found : impossible to test operator against an external reference")
	endif()
endif()
