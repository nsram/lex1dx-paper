# Example consumer of the library: pulls in ../../cpp with add_subdirectory and links lex::lex.
#   cmake -S examples/cpp -B build-examples -DCMAKE_BUILD_TYPE=Release
#   cmake --build build-examples --config Release
#   build-examples/quickstart and build-examples/twelve   (Windows: build-examples\Release\*.exe)
cmake_minimum_required(VERSION 3.21)
project(lex_examples LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()
add_subdirectory(../../cpp lex1d)
add_executable(quickstart quickstart.cpp)
target_link_libraries(quickstart PRIVATE lex::lex)
lex_strict_fp(quickstart)
add_executable(twelve twelve.cpp)
target_link_libraries(twelve PRIVATE lex::lex)
lex_strict_fp(twelve)
