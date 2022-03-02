add_executable(hello_world src/hello_world)
target_link_libraries(hello_world archgenlib Sollya::Sollya)

add_executable(test_expr src/test_expr)
target_link_libraries(test_expr archgenlib Sollya::Sollya HINT_BITINT_BACKEND)

add_custom_command(
    OUTPUT specialization_header.hpp
    COMMAND ./test_expr
    DEPENDS test_expr
    COMMENT "Generating specialization file"
)

add_custom_target(specialization_header DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/specialization_header.hpp )

add_executable(test_expr_2 src/test_expr)
target_link_libraries(test_expr_2 archgenlib HINT_BITINT_BACKEND)
target_include_directories(test_expr_2 PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
target_compile_definitions(test_expr_2 PRIVATE INCLUDE_GENERATED_HEADER)
add_dependencies(test_expr_2 specialization_header)
