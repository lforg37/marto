find_package(Sollya)

add_library(archgenlib SHARED src/sollya_handler)
target_include_directories(archgenlib PUBLIC ${CMAKE_CURRENT_LIST_DIR}/include)
target_link_libraries(archgenlib PRIVATE Sollya::Sollya)

add_subdirectory(examples)
