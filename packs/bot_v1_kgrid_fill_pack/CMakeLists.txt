cmake_minimum_required(VERSION 3.21)
project(kgrid_fill_bench LANGUAGES CXX)
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

# Pristine upstream Grønlund, same flags as scale_cell (not run by this pack).
add_executable(pristine_cell bench/scale_cell.cpp
  vendor/kmeans1d/interval_sum.cpp
  vendor/kmeans1d/kmeans_wilber.cpp)
target_compile_definitions(pristine_cell PRIVATE LEX_PRISTINE)
target_include_directories(pristine_cell PRIVATE vendor/kmeans1d)
lex_strict_fp(pristine_cell)

add_executable(f32_to_f64 bench/f32_to_f64.cpp)
add_executable(txt_to_f64 bench/txt_to_f64.cpp)
add_executable(sort_f64 bench/sort_f64.cpp)
add_executable(score_sse bench/score_sse.cpp
  bench/weighted_wilber/interval_sum.cpp
  bench/weighted_wilber/kmeans_wilber.cpp)
target_include_directories(score_sse PRIVATE bench/weighted_wilber)
lex_strict_fp(score_sse)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  foreach(tgt scale_cell pristine_cell f32_to_f64 txt_to_f64 sort_f64 score_sse)
    target_compile_options(${tgt} PRIVATE -Wall -Wextra)
  endforeach()
endif()
