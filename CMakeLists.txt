cmake_minimum_required(VERSION 3.10)
if(${CMAKE_VERSION} VERSION_GREATER 3.10)
    cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

project(lephare VERSION 1.0 LANGUAGES CXX)

message(STATUS "CMake version: ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}")

# ------------------------------------------------------------------------------
# C++ standard
# ------------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ------------------------------------------------------------------------------
# Build layout
# ------------------------------------------------------------------------------
set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}")
set(CMAKE_BUILD_TYPE "RelWithDebInfo")
set(CMAKE_VERBOSE_MAKEFILE ON)

# ------------------------------------------------------------------------------
# Python
# ------------------------------------------------------------------------------
if(DEFINED ENV{PYTHON_VERSION_REQUIRED})
  message(STATUS "Python version required: $ENV{PYTHON_VERSION_REQUIRED}")
  find_package(
    Python3 $ENV{PYTHON_VERSION_REQUIRED}
    EXACT
    COMPONENTS Development.Module Interpreter
    REQUIRED
  )
else()
  find_package(Python3 COMPONENTS Development.Module Interpreter REQUIRED)
endif()

message(STATUS "Python version: ${Python3_VERSION}")
message(STATUS "Python executable: ${Python3_EXECUTABLE}")
message(STATUS "Python include dir: ${Python3_INCLUDE_DIR}")

# ------------------------------------------------------------------------------
# Sources
# ------------------------------------------------------------------------------
set(SOURCE_DIR "${CMAKE_SOURCE_DIR}/src/lib")
include_directories(${SOURCE_DIR})

set(SOURCES
  ${SOURCE_DIR}/cosmology.cpp
  ${SOURCE_DIR}/ext.cpp
  ${SOURCE_DIR}/flt.cpp
  ${SOURCE_DIR}/globals.cpp
  ${SOURCE_DIR}/keyword.cpp
  ${SOURCE_DIR}/mag.cpp
  ${SOURCE_DIR}/onesource.cpp
  ${SOURCE_DIR}/oneElLambda.cpp
  ${SOURCE_DIR}/opa.cpp
  ${SOURCE_DIR}/PDF.cpp
  ${SOURCE_DIR}/photoz_lib.cpp
  ${SOURCE_DIR}/SED.cpp
  ${SOURCE_DIR}/SEDLib.cpp
)

# ------------------------------------------------------------------------------
# Base compiler flags (always on)
# ------------------------------------------------------------------------------
add_compile_options(
  -Ofast
  -fno-finite-math-only
  -fsigned-zeros
  -fno-associative-math
  -funroll-loops
  -ftree-vectorize
)

# ------------------------------------------------------------------------------
# CPU tuning: dev vs wheels
# ------------------------------------------------------------------------------

option(LEPHARE_NATIVE_OPT "Enable native CPU optimizations" ON)

# Disable native flags when building wheels
if(DEFINED ENV{CIBUILDWHEEL})
  message(STATUS "CIBUILDWHEEL detected: disabling native CPU optimizations")
  set(LEPHARE_NATIVE_OPT OFF CACHE BOOL "" FORCE)
endif()

if(LEPHARE_NATIVE_OPT)
  message(STATUS "Enabling -march=native")
  add_compile_options("-march=native" "-mtune=native")
endif()

# ------------------------------------------------------------------------------
# AVX2 tuning: dev vs wheels for linux
# ------------------------------------------------------------------------------

if(${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
  option(LEPHARE_ENABLE_AVX2 "Enable AVX2 instructions" ON)

  if(DEFINED ENV{CIBUILDWHEEL})
    set(LEPHARE_ENABLE_AVX2 OFF CACHE BOOL "" FORCE)
    add_compile_options("-march=x86-64-v2")
  endif()

  if(LEPHARE_ENABLE_AVX2)
    add_compile_options("-mavx2" "-mno-vzeroupper")
  endif()
endif()


# ------------------------------------------------------------------------------
# Safety guard: forbid native CPU flags and avx2 when building wheels
# ------------------------------------------------------------------------------
if(DEFINED ENV{CIBUILDWHEEL} AND LEPHARE_NATIVE_OPT)
  message(FATAL_ERROR
    "LEPHARE_NATIVE_OPT must be OFF when building wheels (CIBUILDWHEEL detected)"
  )
endif()

if(DEFINED ENV{CIBUILDWHEEL} AND LEPHARE_ENABLE_AVX2)
  message(FATAL_ERROR
    "AVX2 must be disabled when building wheels (CIBUILDWHEEL detected)"
  )
endif()

# ------------------------------------------------------------------------------
# Sanitizers and coverage (env-controlled, unchanged behavior)
# ------------------------------------------------------------------------------
if(DEFINED ENV{SANITIZE_THREAD})
  message(STATUS "Building with thread sanitizer")
  add_compile_options(-g -O2 -fsanitize=thread)
  add_link_options(-fsanitize=thread)
endif()

if(DEFINED ENV{INCLUDE_COVERAGE})
  message(STATUS "Building with coverage instrumentation")
  add_compile_options(-g -O0 -coverage -fprofile-arcs -ftest-coverage)
  add_link_options(-coverage -lgcov -lpthread)
endif()

if(DEFINED ENV{SANITIZE_ADDRESS})
  message(STATUS "Building with address sanitizer")
  add_compile_options(-g -fsanitize=address)
  add_link_options(-fsanitize=address)
endif()

if(DEFINED ENV{INCLUDE_COVERAGE} AND DEFINED ENV{SANITIZE_ADDRESS})
  message(FATAL_ERROR "INCLUDE_COVERAGE and SANITIZE_ADDRESS are incompatible")
endif()

# ------------------------------------------------------------------------------
# Python extension module
# ------------------------------------------------------------------------------
add_subdirectory(extern/pybind11)

pybind11_add_module(
  _lephare
  ${SOURCES}
  ${SOURCE_DIR}/_bindings.cc
)

if(OpenMP_CXX_FOUND)
  target_link_libraries(_lephare PUBLIC OpenMP::OpenMP_CXX)
endif()
