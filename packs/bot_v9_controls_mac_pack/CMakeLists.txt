cmake_minimum_required(VERSION 3.21)
project(v9_controls_bench LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()
set(LEX_BUILD_TESTS OFF CACHE BOOL "Do not build lex tests in this pack" FORCE)
set(LEX_BUILD_CLI OFF CACHE BOOL "Do not build lex_cli in this pack" FORCE)
add_subdirectory(lex1d/cpp)
add_executable(scale_cell bench/scale_cell.cpp
  bench/weighted_wilber/interval_sum.cpp
  bench/weighted_wilber/kmeans_wilber.cpp)
target_include_directories(scale_cell PRIVATE lex1d/cpp/src bench bench/weighted_wilber)
target_link_libraries(scale_cell PRIVATE lex)
lex_strict_fp(scale_cell)
# Pristine upstream Grønlund (vendor/kmeans1d @ 585386d, unmodified) in its own
# binary, same driver and same FP flags as scale_cell. No symbol clash possible.
add_executable(pristine_cell bench/scale_cell.cpp
  vendor/kmeans1d/interval_sum.cpp
  vendor/kmeans1d/kmeans_wilber.cpp)
target_compile_definitions(pristine_cell PRIVATE LEX_PRISTINE)
target_include_directories(pristine_cell PRIVATE vendor/kmeans1d)
lex_strict_fp(pristine_cell)
add_executable(f32_to_f64 bench/f32_to_f64.cpp)
add_executable(txt_to_f64 bench/txt_to_f64.cpp)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(scale_cell PRIVATE -Wall -Wextra)
  target_compile_options(pristine_cell PRIVATE -Wall -Wextra)
  target_compile_options(f32_to_f64 PRIVATE -Wall -Wextra)
  target_compile_options(txt_to_f64 PRIVATE -Wall -Wextra)
endif()
