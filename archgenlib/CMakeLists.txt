find_package(Sollya)

add_library(archgenlib SHARED
    src/expression_rt_tree.cpp
    src/operators.cpp
    src/output_formatter
    src/sollya_handler)
target_include_directories(archgenlib PUBLIC ${CMAKE_CURRENT_LIST_DIR}/include)
target_link_libraries(archgenlib PUBLIC hint)
target_link_libraries(archgenlib PRIVATE Sollya::Sollya)

add_subdirectory(examples)
