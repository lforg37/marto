cmake_minimum_required(VERSION 3.12)
project(PositQuire CXX)

set(VIVADO_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../include)

add_library(libpositquire INTERFACE) 
target_include_directories(libpositquire INTERFACE ${VIVADO_INCLUDE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/includes)

OPTION(BUILD_UNIT_TEST "Build the operators unit tests") 
if(BUILD_UNIT_TEST)
	find_package(Boost 1.55 REQUIRED COMPONENTS unit_test_framework)
	enable_testing()

	add_executable(testConvert_exe tests/testConvert.cpp)
	find_library(SOFTPOSIT_LIB softposit)
	find_path(SOFTPOSIT_H softposit.h)
	target_link_libraries(testConvert_exe libpositquire ${Boost_LIBRARIES}) 
	if(SOFTPOSIT_LIB AND SOFTPOSIT_H)
		MESSAGE(STATUS "softposit found : ${SOFTPOSIT_LIB}, ${SOFTPOSIT_H}")
		target_compile_definitions(testConvert_exe PRIVATE SOFTPOSIT)
		target_include_directories(testConvert_exe PRIVATE ${SOFTPOSIT_H})
		target_link_libraries(testConvert_exe ${SOFTPOSIT_LIB})
	else()
		Message(WARNING "Tests on posit value cannot be run")
	endif()

	add_test(testConvert testConvert_exe)

	#add_executable(IntConstMultShiftAddCostFunction_exe tests/ConstMult/testNodeCost.cpp ${CMAKE_CURRENT_BINARY_DIR}/Factories.cpp ${CMAKE_CURRENT_BINARY_DIR}/VHDLLexer.cpp	)
	#target_include_directories(IntConstMultShiftAddCostFunction_exe PUBLIC ${Boost_INCLUDE_DIR} ${PAGSUITE_INCLUDE_DIR})
	#target_link_libraries(IntConstMultShiftAddCostFunction_exe FloPoCoLib ${Boost_LIBRARIES} ${RPAG_LIB})

	#add_test(IntConstMultShiftAddCost IntConstMultShiftAddCostFunction_exe)
endif()
enable_testing()
