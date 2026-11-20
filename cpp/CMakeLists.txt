cmake_minimum_required(VERSION 3.16)
project(flagalg CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

find_package(OpenMP REQUIRED)
find_path(NAUTY_INCLUDE_DIR nauty.h PATH_SUFFIXES nauty REQUIRED)
find_library(NAUTY_LIBRARY nauty REQUIRED)
find_path(GMP_INCLUDE_DIR gmpxx.h REQUIRED)
find_library(GMP_LIBRARY gmp REQUIRED)
find_library(GMPXX_LIBRARY gmpxx REQUIRED)

add_library(flagalg src/flag.cpp src/theory.cpp src/sdp.cpp src/spec.cpp)
target_include_directories(flagalg PUBLIC include ${NAUTY_INCLUDE_DIR} ${GMP_INCLUDE_DIR})
target_link_libraries(flagalg PUBLIC ${NAUTY_LIBRARY} ${GMPXX_LIBRARY} ${GMP_LIBRARY} OpenMP::OpenMP_CXX)
target_compile_options(flagalg PRIVATE -Wall -Wextra)

add_executable(flagsdp tools/flagsdp.cpp)
target_link_libraries(flagsdp PRIVATE flagalg)

enable_testing()
add_executable(unit_tests tests/unit_tests.cpp)
target_link_libraries(unit_tests PRIVATE flagalg)
add_test(NAME unit_tests COMMAND unit_tests)
