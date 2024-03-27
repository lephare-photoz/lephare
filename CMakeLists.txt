Cmake_minimum_required(VERSION 3.10)

project(lephare VERSION 1.0)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CPP_COMPILER "g++")
set(CMAKE_C_COMPILER "gcc")

message(STATUS "CMAKE_CXX_COMPILER: ${CMAKE_CXX_COMPILER}")
message(STATUS "CMAKE_C_COMPILER: ${CMAKE_C_COMPILER}")

set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}")
set(CMAKE_BINARY_DIR "${CMAKE_SOURCE_DIR}/bin")
set(CMAKE_BUILD_TYPE "RelWithDebInfo")
set(CMAKE_VERBOSE_MAKEFILE ON)

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

# Remove flags that break mac compilation with g++-13
add_compile_options("-Ofast" "-ffast-math" "-fno-finite-math-only" "-fsigned-zeros" "-fno-associative-math" "-fno-inline" "-march=native" "-mtune=native" "-funroll-loops" "-fomit-frame-pointer" "-finline" "-ftree-vectorize")
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
#    add_compile_options()
else()
    add_compile_options("-mno-vzeroupper" "-mavx2" "-floop-unroll-and-jam" )
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
