# Python extension: nanobind bindings over the C++ core in ../cpp (built from the same sources,
# with the same strict floating-point flags (see ../cpp/cmake/LexFloatingPoint.cmake), so results are bit-identical to the C++ library).
cmake_minimum_required(VERSION 3.21)
project(lex_python LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(LEX_CPP_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../cpp" CACHE PATH "C++ core source folder")
if(NOT EXISTS "${LEX_CPP_DIR}/CMakeLists.txt")
  message(FATAL_ERROR "C++ core not found at ${LEX_CPP_DIR}. Build from a full repository checkout.")
endif()
set(LEX_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(LEX_BUILD_CLI OFF CACHE BOOL "" FORCE)
add_subdirectory("${LEX_CPP_DIR}" lex_core EXCLUDE_FROM_ALL)

find_package(Python 3.10 REQUIRED COMPONENTS Interpreter Development.Module)
find_package(nanobind CONFIG REQUIRED)
nanobind_add_module(_core NB_STATIC src/_core.cpp)
target_link_libraries(_core PRIVATE lex::lex)
include("${LEX_CPP_DIR}/cmake/LexFloatingPoint.cmake")
lex_strict_fp(_core)
install(TARGETS _core LIBRARY DESTINATION lex)
