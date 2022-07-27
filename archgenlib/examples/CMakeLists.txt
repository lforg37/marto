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
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/specialization_header.hpp
    COMMAND ${CMAKE_COMMAND} -E env ARCHGENLIB_SPECIALIZATION_HEADER_PATH=${CMAKE_CURRENT_BINARY_DIR}/specialization_header.hpp ./test_expr
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

function(vitis_compile_with_generator target source)
    add_executable(${target}.generator ${source})
    target_link_libraries(${target}.generator archgenlib Sollya::Sollya HINT_BITINT_BACKEND GMPXX::GMPXX)

    add_custom_command(
        OUTPUT specialization_header_${target}.hpp
        COMMAND ${CMAKE_COMMAND} -E env ARCHGENLIB_SPECIALIZATION_HEADER_PATH=${CMAKE_CURRENT_BINARY_DIR}/specialization_header_${target}.hpp ./${target}.generator
        DEPENDS ${target}.generator
        COMMENT "Generating specialization file for " ${source}
    )

    add_custom_target(specialization_header_${target}.hpp DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/specialization_header_${target}.hpp )

    add_custom_command(
        OUTPUT ${target}
        COMMAND ${CMAKE_CXX_COMPILER} --target=vitis_ip-xilinx --vitis-ip-part=xcvu13p-fhga2104-3-e -DBITINT_BACKEND -DBITINT_BACKEND_SIGNED_W1 -DBOOST_ALL_NO_LIB -DINCLUDE_GENERATED_HEADER=${CMAKE_CURRENT_BINARY_DIR}/specialization_header_${target}.hpp -I/storage/gauthier/marto/build/archgenlib/examples -I/storage/gauthier/marto/archgenlib/include -I/storage/gauthier/marto/build/_deps/hint-src/include -fconstexpr-steps=4294967295 -std=c++20 -o ${target}.zip ${CMAKE_CURRENT_SOURCE_DIR}/${source}
        DEPENDS specialization_header_${target}.hpp ${source}
        COMMENT "build a vitis ip for " ${source}
    )
    add_custom_target(${target}.zip ALL DEPENDS ${target})
endfunction()

function(cpu_compile_with_generator target specializatio_header source)
    add_executable(${target} ${source})
    target_link_libraries(${target} archgenlib HINT_BITINT_BACKEND)
    target_include_directories(${target} PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
    target_compile_definitions(${target} PRIVATE INCLUDE_GENERATED_HEADER=${CMAKE_CURRENT_BINARY_DIR}/${specializatio_header})
    target_compile_options(${target} PRIVATE -fconstexpr-steps=4294967295)
    add_dependencies(${target} ${specializatio_header})
endfunction()

vitis_compile_with_generator(test_vitis src/test_vitis.cpp)

vitis_compile_with_generator(test_wave_gen src/test_wave_gen.cpp)
cpu_compile_with_generator(test_wave_gen_cpu specialization_header_test_wave_gen.hpp src/test_wave_gen.cpp)

add_executable(test_literal src/test_literal.cpp)
target_link_libraries(test_literal archgenlib HINT_BITINT_BACKEND)

add_executable(test_lns src/test_lns.cpp)
target_link_libraries(test_lns archgenlib Sollya::Sollya HINT_BITINT_BACKEND GMPXX::GMPXX)

vitis_compile_with_generator(lns src/test_lns.cpp)

add_custom_command(
    OUTPUT wave_gen_float.zip
    COMMAND ${CMAKE_CXX_COMPILER} --target=vitis_ip-xilinx --vitis-ip-part=xcvu13p-fhga2104-3-e -DBITINT_BACKEND -DBITINT_BACKEND_SIGNED_W1 -DBOOST_ALL_NO_LIB -I/storage/gauthier/marto/build/archgenlib/examples -I/storage/gauthier/marto/archgenlib/include -I/storage/gauthier/marto/build/_deps/hint-src/include -fconstexpr-steps=4294967295 -std=c++20 -o wave_gen_float.zip /storage/gauthier/marto/archgenlib/examples/src/test_wave_gen_float.cpp
    DEPENDS src/test_wave_gen_float.cpp
    COMMENT "build a vitis ip"
)
add_custom_target(test_wave_gen_float ALL DEPENDS wave_gen_float.zip )

add_executable(test_wave_gen_float_cpu src/test_wave_gen_float.cpp)
target_link_libraries(test_wave_gen_float_cpu)

add_executable(test_lns_lf src/test_lns_lf.cpp)
target_link_libraries(test_lns_lf archgenlib Sollya::Sollya HINT_BITINT_BACKEND GMPXX::GMPXX)

vitis_compile_with_generator(lns_lf src/test_lns_lf.cpp)

# vitis_compile_with_generator(wave_gen_low_prec src/test_wave_gen_low_prec.cpp)
# cpu_compile_with_generator(test_wave_gen_low_prec_cpu specialization_header_test_wave_gen_low_prec.hpp src/test_wave_gen_low_prec.cpp)

# vitis_compile_with_generator(wave_gen_high_prec src/test_wave_gen_high_prec.cpp)
# cpu_compile_with_generator(test_wave_gen_high_prec_cpu specialization_header_test_wave_gen_high_prec.hpp src/test_wave_gen_high_prec.cpp)

add_executable(sycl_lns src/sycl_lns.cpp)
target_compile_options(sycl_lns PRIVATE -fsycl)
target_link_libraries(sycl_lns archgenlib Sollya::Sollya HINT_BITINT_BACKEND GMPXX::GMPXX)

add_custom_command(
    OUTPUT specialization_header_sycl_lns.hpp
    COMMAND ${CMAKE_COMMAND} -E env ARCHGENLIB_SPECIALIZATION_HEADER_PATH=${CMAKE_CURRENT_BINARY_DIR}/specialization_header_sycl_lns.hpp ./sycl_lns
    DEPENDS sycl_lns
    COMMENT "Generating specialization file"
)

add_custom_target(specialization_header_sycl_lns DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/specialization_header_sycl_lns.hpp )

add_custom_command(
    OUTPUT sycl_app
    COMMAND ${CMAKE_CXX_COMPILER} -fsycl -DBITINT_BACKEND -DBITINT_BACKEND_SIGNED_W1 -DBOOST_ALL_NO_LIB -DINCLUDE_GENERATED_HEADER=${CMAKE_CURRENT_BINARY_DIR}/specialization_header_sycl_lns.hpp -I/storage/gauthier/marto/build/archgenlib/examples -I/storage/gauthier/marto/archgenlib/include -I/storage/gauthier/marto/build/_deps/hint-src/include -fconstexpr-steps=4294967295 -std=c++20 -o sycl_app /storage/gauthier/marto/archgenlib/examples/src/sycl_lns.cpp
    DEPENDS specialization_header_sycl_lns
    COMMENT "build a vitis ip"
)
add_custom_target(sycl_app ALL DEPENDS sycl_app )

vitis_compile_with_generator(wave_gen_6bits src/test_wave_gen_6bits.cpp)
vitis_compile_with_generator(wave_gen_8bits src/test_wave_gen_8bits.cpp)
vitis_compile_with_generator(wave_gen_16bits src/test_wave_gen_16bits.cpp)
vitis_compile_with_generator(wave_gen_22bits src/test_wave_gen_22bits.cpp)
vitis_compile_with_generator(mul_pi src/test_mul_pi.cpp)
