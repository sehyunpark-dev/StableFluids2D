cmake_minimum_required(VERSION 3.10)

project(StableFluids2D VERSION 1.0
        DESCRIPTION "A 2D stable fluids simulation using C++ and OpenGL"
        LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(PROJECT_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include)
set(PROJECT_SOURCE_DIR ${CMAKE_SOURCE_DIR}/src)
set(PROJECT_LIB_DIR ${CMAKE_SOURCE_DIR}/lib)

add_executable(StableFluids2D
    ${PROJECT_SOURCE_DIR}/main.cpp
    ${PROJECT_SOURCE_DIR}/MACGrid2D.cpp
    ${PROJECT_SOURCE_DIR}/Renderer2D.cpp
    ${PROJECT_SOURCE_DIR}/utils.cpp
    ${PROJECT_SOURCE_DIR}/Solver.cpp)

target_include_directories(StableFluids2D PRIVATE 
    ${PROJECT_INCLUDE_DIR}
    ${PROJECT_SOURCE_DIR})

target_link_libraries(StableFluids2D PRIVATE
    glfw3
    glew32s
    opengl32
    gdi32
    user32
    kernel32
    shell32)

target_link_directories(StableFluids2D PRIVATE ${PROJECT_LIB_DIR})

target_compile_definitions(StableFluids2D PRIVATE GLEW_STATIC)

find_package(OpenGL REQUIRED)

find_package(OpenMP)
if(OpenMP_FOUND)
    message(STATUS "OpenMP found, enabling OpenMP support")
    target_compile_options(StableFluids2D PRIVATE ${OpenMP_CXX_FLAGS})
    target_link_libraries(StableFluids2D PRIVATE OpenMP::OpenMP_CXX)
else()
    message(STATUS "OpenMP not found, OpenMP support will not be enabled")
endif()