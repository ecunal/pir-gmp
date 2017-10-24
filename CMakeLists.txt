cmake_minimum_required(VERSION 3.3)
project(pirgmp)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_BUILD_TYPE Release)

find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

include_directories(${PROJECT_SOURCE_DIR})


# GMP
# sudo apt-get install libgmp-dev

#GMP
find_path(GMP_INCLUDES
NAMES
  gmp.h
PATHS
  $ENV{GMPDIR}
  ${INCLUDE_INSTALL_DIR}
)

find_library(GMP_LIBRARIES gmp PATHS $ENV{GMPDIR} ${LIB_INSTALL_DIR})

include_directories(${GMP_INCLUDES})

#OPENMP
FIND_PACKAGE(OpenMP)
IF(OPENMP_FOUND)
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
ENDIF()

# OPTIMIZATION FLAGS
if((CMAKE_CXX_COMPILER_ID MATCHES GNU) OR (CMAKE_CXX_COMPILER_ID MATCHES Clang))
  set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} -std=c++11")
  set(CMAKE_CXX_FLAGS_DEBUG   "-O0 -g3")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
endif()

# GOOGLE BENCHMARK
# git clone https://github.com/google/benchmark.git
# mkdir build && cd build
# cmake .. -DCMAKE_BUILD_TYPE=Release -DBENCHMARK_ENABLE_LTO=true
# make && make install

find_package(benchmark REQUIRED)

set(HEADER_FILES
  Server.h Client.h DamgardJurik.h
)

set(SOURCE_FILES
  Timing.cpp Server.cpp Client.cpp DamgardJurik.cpp
)

add_executable(main
  ${HEADER_FILES}
  ${SOURCE_FILES}
)

add_executable(bench
  bench.cpp
)

target_link_libraries(main ${GTEST_BOTH_LIBRARIES} pthread ${GMP_LIBRARIES})
target_link_libraries(bench ${GTEST_BOTH_LIBRARIES} benchmark pthread ${GMP_LIBRARIES})


