add_executable(hello_world src/hello_world)
target_link_libraries(hello_world archgenlib Sollya::Sollya)

add_executable(test_expr src/test_expr)
target_link_libraries(test_expr archgenlib Sollya::Sollya HINT_BITINT_BACKEND GMPXX::GMPXX)

add_executable(test_func_table src/test_func_table)
target_link_libraries(test_func_table archgenlib Sollya::Sollya HINT_BITINT_BACKEND GMPXX::GMPXX)

add_executable(test_sin src/test_sin)
target_link_libraries(test_sin archgenlib Sollya::Sollya HINT_BITINT_BACKEND GMPXX::GMPXX)

add_executable(fixedpoint_ops src/fixedpoint_operators)
target_link_libraries(fixedpoint_ops archgenlib HINT_BITINT_BACKEND)

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
target_compile_definitions(test_expr_2 PRIVATE INCLUDE_GENERATED_HEADER=${CMAKE_CURRENT_BINARY_DIR}/specialization_header.hpp)
target_compile_options(test_expr_2 PRIVATE -fconstexpr-steps=2560000)
add_dependencies(test_expr_2 specialization_header)

add_executable(test_vitis src/test_vitis.cpp)
target_link_libraries(test_vitis archgenlib Sollya::Sollya HINT_BITINT_BACKEND GMPXX::GMPXX)

add_custom_command(
    OUTPUT specialization_header_vitis.hpp
    COMMAND ./test_vitis
    DEPENDS test_vitis
    COMMENT "Generating specialization file"
)

add_custom_target(specialization_header_vitis DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/specialization_header_vitis.hpp )

add_custom_command(
    OUTPUT function.zip
    COMMAND /storage/gauthier/sycl/build-Release/bin/clang++ --target=vitis_ip-xilinx --vitis-ip-part=xc7vx330t-ffg1157-1 -DBITINT_BACKEND -DBITINT_BACKEND_SIGNED_W1 -DBOOST_ALL_NO_LIB -DINCLUDE_GENERATED_HEADER=${CMAKE_CURRENT_BINARY_DIR}/specialization_header.hpp -I/storage/gauthier/marto/build/archgenlib/examples -I/storage/gauthier/marto/archgenlib/include -I/storage/gauthier/marto/build/_deps/hint-src/include -fconstexpr-steps=4294967295 -std=c++20 -o function.zip /storage/gauthier/marto/archgenlib/examples/src/test_vitis.cpp
    DEPENDS specialization_header_vitis
    COMMENT "build a vitis ip"
)
add_custom_target(test_vitis_2 DEPENDS function.zip )

add_executable(test_wave_gen src/test_wave_gen.cpp)
target_link_libraries(test_wave_gen archgenlib Sollya::Sollya HINT_BITINT_BACKEND GMPXX::GMPXX)

add_custom_command(
    OUTPUT specialization_header_wave_gen.hpp
    COMMAND ./test_wave_gen
    DEPENDS test_wave_gen
    COMMENT "Generating specialization file"
)

add_custom_target(specialization_header_wave_gen DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/specialization_header_wave_gen.hpp )

add_custom_command(
    OUTPUT lns.zip
    COMMAND /storage/gauthier/sycl/build-Release/bin/clang++ --target=vitis_ip-xilinx --vitis-ip-part=xc7vx330t-ffg1157-1 -DBITINT_BACKEND -DBITINT_BACKEND_SIGNED_W1 -DBOOST_ALL_NO_LIB -DINCLUDE_GENERATED_HEADER=${CMAKE_CURRENT_BINARY_DIR}/specialization_header.hpp -I/storage/gauthier/marto/build/archgenlib/examples -I/storage/gauthier/marto/archgenlib/include -I/storage/gauthier/marto/build/_deps/hint-src/include -fconstexpr-steps=4294967295 -std=c++20 -o lns.zip /storage/gauthier/marto/archgenlib/examples/src/test_wave_gen.cpp
    DEPENDS specialization_header_wave_gen
    COMMENT "build a vitis ip"
)
add_custom_target(test_wave_gen_2 DEPENDS lns.zip )
