include(CTest)

find_package(Sollya)
find_package(GMPXX)
find_package(MPFR)

add_library(archgenlib SHARED
    src/expression_tree
    src/fixfunction_multipartite
    src/mpfr_handler
    src/operators
    src/output_formatter
    src/poly_approx
    src/sollya_fix_function
    src/sollya_handler
    src/sollya_operation
    src/utility
    )
target_include_directories(archgenlib PUBLIC ${CMAKE_CURRENT_LIST_DIR}/include)
target_link_libraries(archgenlib PUBLIC hint)
target_link_libraries(archgenlib PRIVATE Sollya::Sollya GMPXX::GMPXX MPFR::MPFR)

if (BUILD_TESTING)
add_subdirectory(examples)
endif()