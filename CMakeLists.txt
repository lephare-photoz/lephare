cmake_minimum_required(VERSION 3.10)
if(${CMAKE_VERSION} VERSION_GREATER 3.10)
    cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()
message(STATUS "CMake version: ${CMAKE_MAJOR_VERSION}")
message(STATUS "CMake minor: ${CMAKE_MINOR_VERSION}")
project(lephare VERSION 1.0)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}")
set(CMAKE_BINARY_DIR "${CMAKE_SOURCE_DIR}/bin")
set(CMAKE_BUILD_TYPE "RelWithDebInfo")
set(CMAKE_VERBOSE_MAKEFILE ON)

if (DEFINED ENV{PYTHON_VERSION_REQUIRED})
  message(STATUS "Python version required: $ENV{PYTHON_VERSION_REQUIRED}")
  find_package(Python3 $ENV{PYTHON_VERSION_REQUIRED} EXACT COMPONENTS Development.Module Interpreter REQUIRED)
else()
  find_package(Python3 COMPONENTS Development.Module Interpreter REQUIRED)
endif()

message(STATUS "Python version: ${Python3_VERSION}")
message(STATUS "Python libraries: ${Python3_LIBRARY}")
message(STATUS "Python executable: ${Python3_EXECUTABLE}")
message(STATUS "Python included_dir: ${Python3_INCLUDE_DIR}")


# Set source directory
set(SOURCE_DIR "./src/lib")
# Tell CMake that headers are also in SOURCE_DIR
include_directories(${SOURCE_DIR})
#provide exhaustive list of source files, as recommended by
#CMake developers
set(SOURCES
  "${SOURCE_DIR}/cosmology.cpp"
  "${SOURCE_DIR}/ext.cpp"
  "${SOURCE_DIR}/flt.cpp"
  "${SOURCE_DIR}/globals.cpp"
  "${SOURCE_DIR}/keyword.cpp"
  "${SOURCE_DIR}/mag.cpp"
  "${SOURCE_DIR}/onesource.cpp"
  "${SOURCE_DIR}/oneElLambda.cpp"
  "${SOURCE_DIR}/opa.cpp"
  "${SOURCE_DIR}/PDF.cpp"
  "${SOURCE_DIR}/photoz_lib.cpp"
  "${SOURCE_DIR}/SED.cpp"
  "${SOURCE_DIR}/SEDLib.cpp"
)

#! May need to remove the -march=native, this may cause problems in the compiled, pip installable, version.
#-0fast enables finite-math-only, but the code is not guaranteed to catch all Nan or inf occurrences, so we need to disable it with -fno-finite-math-only
#-fsigned-zeros is the default in gcc, which does not disable it with -0fast; on the contrary, clang disables it with -0fast so we need to manually enable it.
#tree-vectorize is needed for vectorization, but unroll-loops is for profiling purpose
add_compile_options("-Ofast" "-fno-finite-math-only" "-fsigned-zeros" "-fno-associative-math" "-march=native" "-mtune=native" "-funroll-loops" "-ftree-vectorize")
if(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  #  These break mac compilation with g++-13
  # no-vzeroupper instruct compiler to not optimize control flow from AVX2 to SSE, this is for profiling purpose.
  add_compile_options("-mno-vzeroupper" "-mavx2" )
endif()

# only set these flags if testing presence of race conditions :
# https://github.com/google/sanitizers/wiki/ThreadSanitizerCppManual
# on linux/gcc it seems that one needs to run with LD_PRELOAD=/lib/x86_64-linux-gnu/libtsan.so.0 prepended
# and do a `sudo sysctl vm.mmap_rnd_bits=28` prior to running
if(DEFINED ENV{SANITIZE_THREAD})
  message(STATUS "Building with thread sanitizing enabled.")
  SET(GCC_COVERAGE_COMPILE_FLAGS "-g -O2 -fsanitize=thread")
  SET(GCC_COVERAGE_LINK_FLAGS    "-fsanitize=thread")
  SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")
  SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${GCC_COVERAGE_LINK_FLAGS}")
endif()

if(DEFINED ENV{INCLUDE_COVERAGE} AND DEFINED ENV{SANITIZE_ADDRESS})
  message(FATAL_ERROR "Conflicting build options: INCLUDE_COVERAGE and SANITIZE_ADDRESS")
endif()

# only set these flags if in `INCLUDE_COVERAGE` env var is set (to any value)
if(DEFINED ENV{INCLUDE_COVERAGE})
  message(STATUS "Building with code coverage enabled.")
  SET(GCC_COVERAGE_COMPILE_FLAGS "-g -O0 -coverage -fprofile-arcs -ftest-coverage")
  SET(GCC_COVERAGE_LINK_FLAGS    "-coverage -lgcov -lpthread")
  SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")
  SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${GCC_COVERAGE_LINK_FLAGS}")
endif()

# Enable the address sanitizer with an environment variable.
# This should be on when testing for memory access errors.
# It requires a compiler built with the address sanitizer enabled.
# When running, the library must be preloaded with `LD_PRELOAD=libasan.so`.
# https://github.com/google/sanitizers/wiki/AddressSanitizer
if(DEFINED ENV{SANITIZE_ADDRESS})
  message(STATUS "Building with address sanitizer enabled.")
  SET(GCC_COVERAGE_COMPILE_FLAGS "-g -fsanitize=address")
  SET(GCC_COVERAGE_LINK_FLAGS    "-fsanitize=address")
  SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")
  SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${GCC_COVERAGE_LINK_FLAGS}")
endif()

add_library(lepharelib STATIC ${SOURCES})
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(lepharelib PUBLIC OpenMP::OpenMP_CXX)
endif()

# Build the python module
add_subdirectory(extern/pybind11/)
pybind11_add_module(_lephare ${SOURCES} "${SOURCE_DIR}/_bindings.cc")
# target_compile_definitions(_lephare
#                            PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})

#find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(_lephare PUBLIC OpenMP::OpenMP_CXX)
endif()

