cmake_minimum_required (VERSION 2.8)
project(HOM3)

# Module path
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
set(COMMON_INCLUDES ${PROJECT_SOURCE_DIR}/src)
set(EXT_PROJECTS_DIR ${PROJECT_SOURCE_DIR}/ext)

# Environment
set(CMAKE_OSX_SYSROOT $ENV{SDKROOT})

# Configure compiler:
message("Build type: ${CMAKE_BUILD_TYPE}")
message("The compiler is: ${CMAKE_CXX_COMPILER}")
if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")

  # libcxx:
  find_package(LIBCXX REQUIRED)
  include_directories(SYSTEM ${LIBCXX_INCLUDE_DIR})
  message("libcxx_include: ${LIBCXX_INCLUDE_DIR} | libcxx_lib: ${LIBCXX_LIBRARY}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L${LIBCXX_LIBRARY}")

  # compiler flags:
  set(CMAKE_CXX_FLAGS_LIST
    "-arch x86_64" "-Wall" "-Wextra" "-std=c++1y" "-stdlib=libc++"
    "-I${LIBCXX_INCLUDE_DIR}" "-pedantic" "-Wshadow" "-Woverloaded-virtual"
    "-pedantic-errors" "-Wcast-align" "-Wcomment" "-Wcast-qual"
    "-Wchar-subscripts" "-Wdisabled-optimization" "-Wfloat-equal" "-Wformat=2"
    "-Winvalid-pch" "-Wformat-nonliteral" "-Wformat-security" "-Wformat-y2k"
    "-Wimport" "-Winit-self" "-Winline" "-Wreturn-type" "-Wmissing-braces"
    "-Wmissing-field-initializers" "-Wmissing-include-dirs" "-Wredundant-decls"
    "-Wpacked" "-Wparentheses" "-Wpointer-arith" "-Wsequence-point"
    "-Wsign-compare" "-Wstack-protector" "-Wstrict-aliasing=2" "-Wswitch"
    "-Wswitch-default" "-Wtrigraphs" "-Wuninitialized" "-Wunknown-pragmas"
    "-Wunreachable-code" "-Wunused" "-Wunused-function" "-Wunused-label"
    "-Wunused-parameter" "-Wunused-value" "-Wunused-variable"
    "-Wvariadic-macros" "-Wvolatile-register-var" "-Wwrite-strings"
    "-Woverloaded-virtual" "-Wsign-promo" "-Wstrict-overflow=5"
    "-Wswitch-default" "-DGTEST_USE_OWN_TR1_TUPLE=1"
    "-fdiagnostics-show-template-tree" "-ftemplate-backtrace-limit=0"
    "-Wno-attributes"
    )
  set(CMAKE_CXX_FLAGS_DEBUG_LIST
    "-O0" "-g3" "-fno-inline" "-fstack-protector-all" "-D_FORTIFY_SOURCE=2"
    ) # -D_LIBCPP_DEBUG2=1
  set(CMAKE_CXX_FLAGS_RELEASE_LIST
    "-O3" "-DNDEBUG" "-march=native" "-mtune=native" "-fstrict-aliasing"
    "-fomit-frame-pointer" "-pipe" "-fdata-sections" "-ffunction-sections"
    "-fvectorize" "-fslp-vectorize-aggressive" "-DEIGEN_FAST_MATH"
    "-DEIGEN_NO_DEBUG" "-ffast-math"
    ) # -fno-rtti -fno-exceptions 

  # set flags:
  string(REPLACE ";" " " CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_LIST}")
  string(REPLACE ";" " " CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG_LIST}")
  string(REPLACE ";" " " CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE_LIST}")

  # output flags:
  message("...clang flags set: " ${CMAKE_CXX_FLAGS})
  if("${CMAKE_BUILD_TYPE}" MATCHES "Debug")
    message("...clang optimization flags set: ${CMAKE_CXX_FLAGS_DEBUG}")
  elseif("${CMAKE_BUILD_TYPE}" MATCHES "Release")
    message("...clang optimization flags set: ${CMAKE_CXX_FLAGS_RELEASE}")
  else()
    message(FATAL_ERROR "Unknown build type. Options are \"Debug\" and \"Release\". Quit.")
  endif()
  message("...linker flags set: ${CMAKE_EXE_LINKER_FLAGS}")
else()
  message(FATAL_ERROR "Unknown compiler, quit.")
endif()

# Configure valgrind
set(CMAKE_MEMORYCHECK_COMMAND "valgrind")
# set(CMAKE_MEMORYCHECK_COMMAND_OPTIONS "--trace-children=yes --quiet
# --tool=memcheck --leak-check=yes --show-reachable=yes --num-callers=100
# --verbose --demangle=yes")
set(CMAKE_MEMORYCHECK_COMMAND_OPTIONS "--quiet --tool=memcheck --error-exitcode=1")
function(add_memcheck_test name binary)
  set(memcheck_command "${CMAKE_MEMORYCHECK_COMMAND} ${CMAKE_MEMORYCHECK_COMMAND_OPTIONS}")
  separate_arguments(memcheck_command)
  add_test(${name} ${binary} ${ARGN})
  add_test(memcheck_${name} ${memcheck_command} ./${binary} ${ARGN})
endfunction(add_memcheck_test)

# Eigen
find_package(Eigen3 REQUIRED)
set(CMAKE_INCLUDE_SYSTEM_FLAG_CXX "-isystem ")
include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})

# gTest (and enables unit testing)
add_subdirectory(${EXT_PROJECTS_DIR}/gtest)
include_directories(${GTEST_INCLUDE_DIRS} ${COMMON_INCLUDES})
enable_testing(true)

set(HOM3_BASE_LIBS ${Boost_LIBRARIES} ${MPI_LIBRARIES} ${PNETCDF_LIBRARY})
set(HOM3_TESTING_LIBS ${GTEST_LIBS_DIR}/libgtest.a ${GTEST_LIBS_DIR}/libgtest_main.a)
set(HOM3_TESTING_INCLUDES ${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
function(add_hom3_test name)
  include_directories(${HOM3_TESTING_INCLUDES})
  add_executable(${name} ${name}_test.cpp)
  add_dependencies(${name} googletest)
  target_link_libraries(${name} ${HOM3_TESTING_LIBS})
  add_test(${name} ${name})
endfunction(add_hom3_test)

# MPI
find_package(MPI REQUIRED)
set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
include_directories(MPI_INCLUDE_PATH)

# boost
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME ON)
find_package(Boost 1.53 COMPONENTS mpi serialization filesystem system REQUIRED)
include_directories(SYSTEM ${Boost_INCLUDE_DIRS})

# Doxygen (optional) - Enables "make docs"
option(BUILD_DOCUMENTATION "Use Doxygen to create the HTML documentation" ON)
if(BUILD_DOCUMENTATION)
  find_package(Doxygen)
  if (NOT DOXYGEN_FOUND)
    message(FATAL_ERROR "Doxygen can't be found.")
  endif()
  configure_file(Doxyfile.in ${PROJECT_BINARY_DIR}/Doxyfile @ONLY IMMEDIATE)
  add_custom_target (docs
    ${DOXYGEN_EXECUTABLE}
    ${PROJECT_BINARY_DIR}/Doxyfile
    SOURCES
    ${PROJECT_BINARY_DIR}/Doxyfile)
endif()

# Subdirectories:
add_subdirectory (./src/misc/tests)
add_subdirectory (./src/containers/sequential/tests)
add_subdirectory (./src/containers/hierarchical/tests)
add_subdirectory (./src/grid/)
add_subdirectory (./src/output/)
add_subdirectory (./src/solver/fv)
add_subdirectory (./src/geometry/)
