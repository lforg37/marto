cmake_minimum_required(VERSION 3.10)
project(Marto CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(LIBNAME marto)
set(INCLUDE_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/includes)

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

	find_path(MPFR mpfr.h)

	if(MPFR)
		MESSAGE(STATUS "MPFR found : ${MPFR}")
		add_executable(testKulisch_exe tests/testKulisch.cpp)
		target_compile_definitions(testKulisch_exe PRIVATE MPFR)
		target_compile_options(testKulisch_exe PRIVATE -O3 -lmpfr)
		target_link_libraries(testKulisch_exe PRIVATE ${LIBNAME} ${Boost_LIBRARIES} VIVADO_AP_INT -lmpfr)
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
