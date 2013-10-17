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
  set(CMAKE_EXE_LINKER_FLAGS " ${CMAKE_EXE_LINKER_FLAGS} -L${LIBCXX_LIBRARY} -lc++")
  message("linker flags: ${CMAKE_EXE_LINKER_FLAGS}")

  # link using clang or the sanitizers _might_ fail to link:
  set(CMAKE_CXX_LINK_EXECUTABLE=${CMAKE_CXX_COMPILER})

  # compiler flags:
  set(CMAKE_CXX_FLAGS_LIST
    "-arch x86_64" "-Wall" "-Wextra" "-std=c++1y" "-stdlib=libc++"
    "-pedantic" "-Wshadow" "-Woverloaded-virtual"
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
    "-Wno-attributes" "-DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_TBB"
    )
  set(CMAKE_CXX_FLAGS_DEBUG_LIST
    "-O1" "-g3" "-fstack-protector-all" "-D_FORTIFY_SOURCE=2"
    ) # -D_LIBCPP_DEBUG2=1
  set(CMAKE_CXX_FLAGS_RELEASE_LIST
    "-O3" "-DNDEBUG" "-march=native" "-mtune=native" "-fstrict-aliasing"
    "-fomit-frame-pointer" "-pipe" "-fdata-sections" "-ffunction-sections"
    "-fvectorize" "-fslp-vectorize-aggressive" "-DEIGEN_FAST_MATH"
    "-DEIGEN_NO_DEBUG" "-ffast-math"
    ) # -fno-rtti -fno-exceptions
  set(CMAKE_CXX_FLAGS_ASAN_LIST
    ${CMAKE_CXX_FLAGS_RELEASE_LIST}
    "-fsanitize=address" "-fno-omit-frame-pointer"
    )
  set(CMAKE_CXX_FLAGS_MSAN_LIST
    ${CMAKE_CXX_FLAGS_RELEASE_LIST}
    "-fsanitize=memory" "-fno-omit-frame-pointer"
    "-fsanitize-memory-track-origins"
    )
  set(CMAKE_CXX_FLAGS_TSAN_LIST
    ${CMAKE_CXX_FLAGS_RELEASE_LIST}
    "-fsanitize=thread" "-fno-omit-frame-pointer"
    "-fsanitize-memory-track-origins"
    )

  # set flags:
  string(REPLACE ";" " " CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_LIST}")
  string(REPLACE ";" " " CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG_LIST}")
  string(REPLACE ";" " " CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE_LIST}")
  string(REPLACE ";" " " CMAKE_CXX_FLAGS_ASAN "${CMAKE_CXX_FLAGS_ASAN_LIST}")
  string(REPLACE ";" " " CMAKE_CXX_FLAGS_MSAN "${CMAKE_CXX_FLAGS_MSAN_LIST}")
  string(REPLACE ";" " " CMAKE_CXX_FLAGS_TSAN "${CMAKE_CXX_FLAGS_TSAN_LIST}")

  # output flags:
  message("...clang flags set: " ${CMAKE_CXX_FLAGS})
  if("${CMAKE_BUILD_TYPE}" MATCHES "Debug")
    message("...clang optimization flags set: ${CMAKE_CXX_FLAGS_DEBUG}")
  elseif("${CMAKE_BUILD_TYPE}" MATCHES "Release")
    message("...clang optimization flags set: ${CMAKE_CXX_FLAGS_RELEASE}")
  elseif("${CMAKE_BUILD_TYPE}" MATCHES "Asan")
    message("...clang optimization flags set: ${CMAKE_CXX_FLAGS_ASAN}")
    set(CMAKE_EXE_LINKER_FLAGS " ${CMAKE_EXE_LINKER_FLAGS} -L${LIBCXX_LIBRARY} -lc++")
  elseif("${CMAKE_BUILD_TYPE}" MATCHES "Msan")
    message("...clang optimization flags set: ${CMAKE_CXX_FLAGS_MSAN}")
  elseif("${CMAKE_BUILD_TYPE}" MATCHES "Tsan")
    message("...clang optimization flags set: ${CMAKE_CXX_FLAGS_TSAN}")
  else()
    message(FATAL_ERROR "Unknown build type. Options are \"Debug\" and \"Release\". Quit.")
  endif()
  message("...linker flags set: ${CMAKE_EXE_LINKER_FLAGS}")
else()
  message(FATAL_ERROR "Unknown compiler, quit.")
endif()

# Eigen
find_package(Eigen3 REQUIRED)
set(CMAKE_INCLUDE_SYSTEM_FLAG_CXX "-isystem ")
include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})

# gTest (and enables unit testing)
add_subdirectory(${EXT_PROJECTS_DIR}/gtest)
include_directories(${GTEST_INCLUDE_DIRS} ${COMMON_INCLUDES})
set(GTEST_LIBS ${GTEST_LIBS_DIR}/libgtest.a ${GTEST_LIBS_DIR}/libgtest_main.a)
enable_testing(true)

# MPI
find_package(MPI REQUIRED)
set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
include_directories(MPI_INCLUDE_PATH)

# HDF5
find_package(HDF5 REQUIRED)
set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${HDF5_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${HDF5_LINK_FLAGS})
include_directories(HDF5_INCLUDE_PATH)

# HDF5-Lite
find_package(HDF5LITE REQUIRED)
set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${HDF5LITE_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${HDF5LITE_LINK_FLAGS})
include_directories(HDF5LITE_INCLUDE_PATH)

# boost
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME ON)
find_package(Boost 1.53 COMPONENTS mpi serialization filesystem system REQUIRED)
include_directories(SYSTEM ${Boost_INCLUDE_DIRS})

# Thrust
find_package(Thrust REQUIRED)
include_directories(SYSTEM ${THRUST_INCLUDE_DIR})

find_package(TBB REQUIRED)
set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${TBB_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${TBB_LINK_FLAGS})
include_directories(TBB_INCLUDE_PATH)

set(HOM3_BASE_LIBS ${Boost_LIBRARIES} ${MPI_LIBRARIES} ${HDF5_LIBRARY} ${HDF5LITE_LIBRARY} ${TBB_LIBRARIES})
set(HOM3_TESTING_LIBS ${HOM3_BASE_LIBS} ${GTEST_LIBS})
set(HOM3_TESTING_INCLUDES ${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
function(add_hom3_test name)
  include_directories(${HOM3_TESTING_INCLUDES})
  add_executable(${name} ${name}_test.cpp)
  add_dependencies(${name} googletest)
  target_link_libraries(${name} ${HOM3_TESTING_LIBS})
  add_test(${name} ${name})
endfunction(add_hom3_test)

function(add_hom3_mpi_test name no_mpi_proc)
  include_directories(${HOM3_TESTING_INCLUDES})
  add_executable(${name} ${name}_test.cpp)
  add_dependencies(${name} googletest)
  target_link_libraries(${name} ${HOM3_TESTING_LIBS})
  set(test_parameters -np ${no_mpi_proc} "./${name}")
  add_test(NAME ${name} COMMAND "mpirun" ${test_parameters})
endfunction(add_hom3_mpi_test)

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

# Clang-format - enables "make format"
set(FORMAT_SCRIPT ${PROJECT_BINARY_DIR}/tools/style-check.sh)
add_custom_target (style ${FORMAT_SCRIPT})

# Clean-up - enables "make clean"
set(CLEANUP_SCRIPT ${PROJECT_BINARY_DIR}/tools/cleanup.sh)
add_custom_target (clean ${CLEANUP_SCRIPT})

# Subdirectories:
add_subdirectory (./src)
