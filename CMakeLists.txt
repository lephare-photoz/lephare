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
# Remove flags that break mac compilation with g++-13
add_compile_options("-Ofast" "-ffast-math" "-fno-finite-math-only" "-fsigned-zeros" "-fno-associative-math" "-fno-inline" "-march=native" "-mtune=native" "-funroll-loops" "-fomit-frame-pointer" "-finline" "-ftree-vectorize")
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(ENV{CC} "/usr/local/opt/llvm/bin/clang")
    set(ENV{CXX} "/usr/local/opt/llvm/bin/clang++")
    set(ENV{LDFLAGS} "-L/usr/local/opt/llvm/lib")
    set(ENV{CPPFLAGS} "-I/usr/local/opt/llvm/include")
    set(CMAKE_CPP_COMPILER "clang++")
    set(CMAKE_C_COMPILER "clang")
else()
    add_compile_options("-mno-vzeroupper" "-mavx2" "-floop-unroll-and-jam" )
endif()

# only set these flags if in `INCLUDE_COVERAGE` env var is set (to any value)
if(DEFINED ENV{INCLUDE_COVERAGE})
  message(STATUS "Building with code coverage enabled.")
  SET(GCC_COVERAGE_COMPILE_FLAGS "-g -O0 -coverage -fprofile-arcs -ftest-coverage")
  SET(GCC_COVERAGE_LINK_FLAGS    "-coverage -lgcov")
  SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}" )
  SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${GCC_COVERAGE_LINK_FLAGS}" )
endif()

add_library(lepharelib STATIC ${SOURCES})
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(lepharelib PUBLIC OpenMP::OpenMP_CXX)
endif()

# Build the LePhare executables
add_executable(filter
  "${SOURCE_DIR}/filter.cpp"
  )
target_link_libraries(filter lepharelib)

add_executable(sedtolib
  "${SOURCE_DIR}/sedtolib.cpp"  
  )
target_link_libraries(sedtolib lepharelib)

add_executable(mag_gal
  "${SOURCE_DIR}/mag_gal.cpp"  
  )
target_link_libraries(mag_gal lepharelib)

add_executable(zphota
  "${SOURCE_DIR}/zphota.cpp"  
  )
target_link_libraries(zphota lepharelib)

add_executable(filter_extinc
  "${SOURCE_DIR}/filter_extinc.cpp"
  )
target_link_libraries(filter_extinc lepharelib)

# Build the python module
add_subdirectory(extern/pybind11/)
pybind11_add_module(_lephare ${SOURCES} "${SOURCE_DIR}/_bindings.cc")
# target_compile_definitions(_lephare
#                            PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})

#find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(_lephare PUBLIC OpenMP::OpenMP_CXX)
endif()

install(TARGETS filter RUNTIME DESTINATION bin)
install(TARGETS sedtolib RUNTIME DESTINATION bin)
install(TARGETS mag_gal RUNTIME DESTINATION bin)
install(TARGETS zphota RUNTIME DESTINATION bin)
install(TARGETS filter_extinc RUNTIME DESTINATION bin)
