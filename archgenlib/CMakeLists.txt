include(CTest)

find_package(Sollya)
find_package(GMPXX)
find_package(MPFR)

add_library(archgenlib SHARED
    src/expression_tree.cpp
    src/fixfunction_multipartite.cpp
    src/mpfr_handler.cpp
    src/operators.cpp
    src/output_formatter.cpp
    src/poly_approx.cpp
    src/sollya_fix_function.cpp
    src/sollya_handler.cpp
    src/sollya_operation.cpp
    src/utility.cpp
    )
target_include_directories(archgenlib PUBLIC ${CMAKE_CURRENT_LIST_DIR}/include)
target_link_libraries(archgenlib PUBLIC hint)
target_link_libraries(archgenlib PRIVATE Sollya::Sollya GMPXX::GMPXX MPFR::MPFR)

if (BUILD_TESTING)
add_subdirectory(examples)
endif()
