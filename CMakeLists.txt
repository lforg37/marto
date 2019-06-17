cmake_minimum_required(VERSION 3.10)
project(Marto CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(LIBNAME marto)
set(GENERATED_HEADER_LIBNAME generated_header)
set(INCLUDE_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/includes)
set(GENERATED_INCLUDE_ROOT ${CMAKE_CURRENT_BINARY_DIR}/includes)
set(GENERATED_HEADER)
set(GENERATED_LIST ${CMAKE_CURRENT_BINARY_DIR}/list.txt)

set(SUBDIR_LIST_FILE ${CMAKE_CURRENT_SOURCE_DIR}/IncDirList.txt)
file(STRINGS ${SUBDIR_LIST_FILE} SUBDIR_LIST)

find_package(libhint CONFIG REQUIRED)

if(NOT DEFINED AP_INT_INCLUDE_DIR)
	MESSAGE(FATAL_ERROR "Please provide the path to ap_int include directory in AP_INT_INCLUDE_DIR")
endif()

set(AP_INT_INCLUDE_DIR ${AP_INT_INCLUDE_DIR} CACHE PATH "Path to ap_int include directory" FORCE)
macro(add_subdir dirname)
	message(STATUS "Call add_subdir for dirname " ${dirname})
	set(dependency ${GENERATED_INCLUDE_ROOT}/${LIBNAME}/${dirname}.hpp)
	if(
		NOT EXISTS ${dependency} 
		OR 
			${INCLUDE_ROOT}/${dirname}/list.txt 
			IS_NEWER_THAN 
			${dependency}
	)
		execute_process(COMMAND
			${CMAKE_COMMAND} 
			-DPREFIX=${dirname}/
			-DOUTPUT=${dependency}
			-DMERGE_LIST=${INCLUDE_ROOT}/${dirname}/list.txt
			-P ${CMAKE_CURRENT_SOURCE_DIR}/merge_headers.cmake
		)
	endif()
	#List the includes in per-directory basis
	add_custom_command(
		OUTPUT ${dependency}
		COMMAND ${CMAKE_COMMAND} 
			-DPREFIX=${dirname}/
			-DOUTPUT=${dependency}
			-DMERGE_LIST=${INCLUDE_ROOT}/${dirname}/list.txt
			-P ${CMAKE_CURRENT_SOURCE_DIR}/merge_headers.cmake
		MAIN_DEPENDENCY ${INCLUDE_ROOT}/${dirname}/list.txt
	)
	set(GENERATED_HEADER ${GENERATED_HEADER} ${dependency})
	target_sources(${GENERATED_HEADER_LIBNAME} PUBLIC ${dependency})
endmacro(add_subdir)

add_library(VIVADO_AP_INT INTERFACE IMPORTED)
target_include_directories(VIVADO_AP_INT INTERFACE ${AP_INT_INCLUDE_DIR})

add_library(${GENERATED_HEADER_LIBNAME} STATIC)
add_library(${LIBNAME} INTERFACE) 
set_target_properties(${GENERATED_HEADER_LIBNAME}
    PROPERTIES LINKER_LANGUAGE CXX
   )
target_link_libraries(${LIBNAME} INTERFACE hint ${GENERATED_HEADER_LIBNAME})

file(REMOVE ${GENERATED_INCLUDE_ROOT}/${LIBNAME}.hpp)

foreach(header_dirname ${SUBDIR_LIST})
	get_filename_component(dirname ${header_dirname} NAME_WE)
	add_subdir(${dirname})
endforeach()

execute_process(COMMAND
			${CMAKE_COMMAND} 
			-DPREFIX=${LIBNAME}/ 
			-DOUTPUT=${GENERATED_INCLUDE_ROOT}/${LIBNAME}.hpp 
			-DMERGE_LIST=${SUBDIR_LIST_FILE}
			-P ${CMAKE_CURRENT_SOURCE_DIR}/merge_headers.cmake
		)
add_custom_command(
		OUTPUT ${GENERATED_INCLUDE_ROOT}/${LIBNAME}.hpp 
		COMMAND ${CMAKE_COMMAND} 
			-DPREFIX=${LIBNAME}/ 
			-DOUTPUT=${GENERATED_INCLUDE_ROOT}/${LIBNAME}.hpp 
			-DMERGE_LIST=${SUBDIR_LIST_FILE}
			-P ${CMAKE_CURRENT_SOURCE_DIR}/merge_headers.cmake
		MAIN_DEPENDENCY ${INCLUDE_ROOT}/${dirname}/list.txt
	)
set(GENERATED_HEADER ${GENERATED_HEADER} ${GENERATED_INCLUDE_ROOT}/${LIBNAME}.hpp) 
target_sources(${GENERATED_HEADER_LIBNAME} PUBLIC ${GENERATED_INCLUDE_ROOT}/${LIBNAME}.hpp)

target_include_directories(${GENERATED_HEADER_LIBNAME} PUBLIC ${GENERATED_INCLUDE_ROOT})
target_include_directories(${LIBNAME} INTERFACE 
	$<BUILD_INTERFACE:${INCLUDE_ROOT}> 
	$<INSTALL_INTERFACE:include/marto>
)

install(DIRECTORY ${INCLUDE_ROOT}/ DESTINATION include/marto FILES_MATCHING PATTERN "*.hpp" PATTERN "*.ipp") 
install(DIRECTORY ${GENERATED_INCLUDE_ROOT}/ DESTINATION include FILES_MATCHING PATTERN "*.hpp" PATTERN "*.ipp") 

OPTION(BUILD_UNIT_TEST "Build the operators unit tests") 
if(BUILD_UNIT_TEST)
	find_package(Boost 1.55 REQUIRED COMPONENTS unit_test_framework)
	enable_testing()

	set(TESTDATAROOT ${CMAKE_CURRENT_SOURCE_DIR}/testdata)

	add_executable(testConvert_exe tests/testConvert.cpp)
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
		target_compile_options(testAritOps_exe PRIVATE -O3 -fopenmp)
		target_link_libraries(testAritOps_exe ${LIBNAME} ${Boost_LIBRARIES} ${SOFTPOSIT_LIB} VIVADO_AP_INT -fopenmp)
		target_compile_definitions(testAritOps_exe PRIVATE TESTDATA_ROOT="${TESTDATAROOT}" VIVADO_BACKEND)
		add_test(testAritOps testAritOps_exe)
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

	find_path(GMP_INCLUDE_DIR gmpxx.h)
	find_library(GMP gmp)
	find_library(GMPXX gmpxx)

	if (NOT (GMPXX-NOTFOUND OR GMP-NOTFOUND OR GMP_INCLUDE_DIR-NOTFOUND))
		set(GMPTARGET-OK ON)
		message(STATUS "GMP found")
	else()
		set(GMPTARGET-NOTFOUND ON)
		message(STATUS "GMP Not found")
	endif()

	if(GMPTARGET-OK)
		add_executable(testIEEEAdderGMP_exe tests/testIEEEAdderGMP.cpp)
		target_compile_definitions(testIEEEAdderGMP_exe PRIVATE FASTSIM_BACKEND)
		target_link_libraries(testIEEEAdderGMP_exe PRIVATE ${LIBNAME} ${Boost_LIBRARIES} ${GMP} ${GMPXX})
		target_compile_definitions(testIEEEAdderGMP_exe PRIVATE TESTDATA_ROOT="${TESTDATAROOT}")
		add_test(testIEEEAdderGMP testIEEEAdderGMP_exe)
	endif()
endif()
