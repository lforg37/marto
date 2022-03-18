add_executable(hello_world src/hello_world)
target_link_libraries(hello_world archgenlib Sollya::Sollya)

add_executable(test_expr src/test_expr)
target_link_libraries(test_expr archgenlib Sollya::Sollya HINT_BITINT_BACKEND GMPXX::GMPXX)

add_executable(test_func_table src/test_func_table)
target_link_libraries(test_func_table archgenlib Sollya::Sollya HINT_BITINT_BACKEND GMPXX::GMPXX)

add_executable(test_sin src/test_sin)
target_link_libraries(test_sin archgenlib Sollya::Sollya HINT_BITINT_BACKEND GMPXX::GMPXX)

add_executable(test_sin_mp src/test_sin_mp)
target_link_libraries(test_sin_mp archgenlib Sollya::Sollya HINT_BITINT_BACKEND GMPXX::GMPXX)
#target_compile_options(test_sin_mp PRIVATE "-fsanitize=address,undefined")
#target_link_options(test_sin_mp PRIVATE "-fsanitize=address,undefined")
add_test(multipartite_sin test_sin_mp)

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
target_compile_options(test_expr_2 PRIVATE -fconstexpr-steps=2560000)
add_dependencies(test_expr_2 specialization_header)
