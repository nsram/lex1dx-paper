cmake_minimum_required(VERSION 3.16)
project(v11_ckmeans_matched LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

add_executable(sort_f64 bench/sort_f64.cpp)
add_executable(txt_to_f64 bench/txt_to_f64.cpp)
add_executable(score_sse bench/score_sse.cpp
  bench/weighted_wilber/interval_sum.cpp
  bench/weighted_wilber/kmeans_wilber.cpp)
target_include_directories(score_sse PRIVATE bench/weighted_wilber)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  foreach(tgt sort_f64 txt_to_f64 score_sse)
    target_compile_options(${tgt} PRIVATE -Wall -Wextra -ffp-contract=off)
  endforeach()
endif()
