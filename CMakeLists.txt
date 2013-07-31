################################################################################
## This is the main HOM3 configuration file
################################################################################
cmake_minimum_required (VERSION 2.8)
project(HOM3)

# Select the build type
# TODO: this should be done by passing an option to cmake
#set(CMAKE_BUILD_TYPE Release_OMP)
#set(CMAKE_BUILD_TYPE Release)
#set(CMAKE_BUILD_TYPE Debug)
#set(CMAKE_BUILD_TYPE Debug_Extra)

# Module path
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# Environment
set(CMAKE_OSX_SYSROOT $ENV{SDKROOT})

# Compiler flags
MESSAGE("Build type: ${CMAKE_BUILD_TYPE}")
IF("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
  MESSAGE("Using gcc-4.8...")
  set(CMAKE_CXX_FLAGS "-arch x86_64 -Wextra -Wall -std=c++11 -pedantic -Wshadow -Woverloaded-virtual -pedantic-errors -Wcast-align -Wcomment -Wcast-qual -Wchar-subscripts  -Wdisabled-optimization -Wfloat-equal -Wformat -Wformat=2 -Winvalid-pch -Wformat-nonliteral -Wformat-security -Wformat-y2k -Wimport -Winit-self -Winline -Wreturn-type -Wunsafe-loop-optimizations -Wmissing-braces -Wmissing-field-initializers -Wmissing-include-dirs -Wredundant-decls -Wpacked -Wparentheses -Wpointer-arith -Wsequence-point -Wsign-compare -Wstack-protector -Wstrict-aliasing -Wstrict-aliasing=2 -Wswitch -Wswitch-default  -Wtrigraphs -Wuninitialized -Wunknown-pragmas -Wunreachable-code -Wunused -Wunused-function -Wunused-label -Wunused-parameter -Wunused-value -Wunused-variable -Wvariadic-macros -Wvolatile-register-var -Wwrite-strings -Wlogical-op -Wnoexcept -Woverloaded-virtual -Wsign-promo -Wstrict-null-sentinel -Wstrict-overflow=1 -Wswitch-default")
  set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g3 -fno-inline -fstack-protector-all -D_FORTIFY_SOURCE=2 -D_GLIBCXX_DEBUG")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native -mtune=native -fstrict-aliasing -funroll-loops -fno-exceptions -fomit-frame-pointer -pipe -fdata-sections -ffunction-sections  -fprefetch-loop-arrays -funsafe-loop-optimizations -ffast-math  -flto -fwhole-program -floop-interchange -floop-strip-mine -ftree-loop-linear -floop-block -floop-flatten -floop-parallelize-all -mno-avx")
# -fno-rtti -fuse-linker-plugin -Wl,--gc-sections
  set(CMAKE_CXX_FLAGS_RELEASE_OMP "-O3 -DNDEBUG -march=native -mtune=native -fstrict-aliasing -funroll-loops -fno-exceptions -fomit-frame-pointer -pipe -fdata-sections -ffunction-sections -fprefetch-loop-arrays -funsafe-loop-optimizations -ffast-math -flto -fwhole-program -floop-interchange -floop-strip-mine -ftree-loop-linear -floop-block -floop-flatten -floop-parallelize-all -fopenmp -D_GLIBCXX_PARALLEL -mno-avx")
# -fno-rtti -fuse-linker-plugin -Wl,--gc-sections
  set(CMAKE_CXX_FLAGS_DEBUG_EXTRA "-Wpadded -Wundef -Weffc++ -Wabi -Wsign-conversion -Wctor-dtor-privacy -Wconversion -Wmissing-format-attribute -Wold-style-cast -Wswitch-enum -Wmissing-noreturn -fstack-protector-all -D_FORTIFY_SOURCE=2 -g3 -fno-inline -O0")
  MESSAGE("...gcc flags set.")
ELSEIF("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
  # libcxx
  find_package(LIBCXX REQUIRED)
  include_directories(SYSTEM ${LIBCXX_INCLUDE_DIR})
  MESSAGE("Using Clang...")
  MESSAGE("libcxx_include: ${LIBCXX_INCLUDE_DIR} | libcxx_lib: ${LIBCXX_LIBRARY}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L${LIBCXX_LIBRARY}")
  set(CMAKE_CXX_FLAGS "-arch x86_64 -Wall -Wextra -std=c++11 -stdlib=libc++ -I${LIBCXX_INCLUDE_DIR} -pedantic -Wshadow -Woverloaded-virtual -pedantic-errors -Wcast-align -Wcomment -Wcast-qual  -Wchar-subscripts  -Wdisabled-optimization -Wfloat-equal -Wformat -Wformat=2 -Winvalid-pch -Wformat-nonliteral -Wformat-security -Wformat-y2k -Wimport  -Winit-self  -Winline -Wreturn-type  -Wmissing-braces -Wmissing-field-initializers -Wmissing-include-dirs  -Wredundant-decls -Wpacked -Wparentheses  -Wpointer-arith -Wsequence-point  -Wsign-compare  -Wstack-protector -Wstrict-aliasing -Wstrict-aliasing=2 -Wswitch -Wswitch-default  -Wtrigraphs -Wuninitialized -Wunknown-pragmas  -Wunreachable-code -Wunused -Wunused-function -Wunused-label  -Wunused-parameter -Wunused-value  -Wunused-variable -Wvariadic-macros -Wvolatile-register-var  -Wwrite-strings -Woverloaded-virtual -Wsign-promo -Wstrict-overflow=5 -Wswitch-default -DGTEST_USE_OWN_TR1_TUPLE=1 -fdiagnostics-show-template-tree -ftemplate-backtrace-limit=0 -Wno-attributes")
  set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g3 -fno-inline -fstack-protector-all -D_FORTIFY_SOURCE=2 -D_GLIBCXX_DEBUG")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native -mtune=native -fstrict-aliasing -fomit-frame-pointer -pipe -fdata-sections -ffunction-sections -fvectorize -fslp-vectorize-aggressive -DEIGEN_FAST_MATH -DEIGEN_NO_DEBUG") # -fno-rtti -fno-exceptions -ffast-math 
  MESSAGE("...clang flags set: ${CMAKE_CXX_FLAGS}")
  IF("${CMAKE_BUILD_TYPE}" MATCHES "Debug")
    MESSAGE("...clang optimization flags set: ${CMAKE_CXX_FLAGS_DEBUG}")
  ELSEIF("${CMAKE_BUILD_TYPE}" MATCHES "Release")
    MESSAGE("...clang optimization flags set: ${CMAKE_CXX_FLAGS_RELEASE}")
  ELSE()
    MESSAGE(FATAL_ERROR "Unknown build type. Options are \"Debug\" and \"Release\". Quit.")
  ENDIF()
  MESSAGE("...linker flags set: ${CMAKE_EXE_LINKER_FLAGS}")
ELSE()
  MESSAGE(FATAL_ERROR "Unknown compiler. Options are \"g++-4.8\" and \"clang++\". Quit.")
ENDIF()

################################################################################
# Configure third party libraries:

# Eigen
find_package(Eigen3 REQUIRED)
set(CMAKE_INCLUDE_SYSTEM_FLAG_CXX "-isystem ")
include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})

# gTest (and enables unit testing)
add_subdirectory(gtest-1.6.0)
enable_testing(true)

# MPI (this removes the need of using unflexible wrappers like e.g. mpicxx)
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

# pnetcdf
#find_package(PNETCDF REQUIRED)
#include_directories(SYSTEM ${PNETCDF_INCLUDE_DIR})



SET(HOM3_BASE_LIBS ${Boost_LIBRARIES} ${MPI_LIBRARIES} ${PNETCDF_LIBRARY})

# Doxygen (optional) - Enables "make docs"
option(BUILD_DOCUMENTATION "Use Doxygen to create the HTML documentation" ON)
if(BUILD_DOCUMENTATION)
  find_package(Doxygen)
  if (NOT DOXYGEN_FOUND)
    message(FATAL_ERROR "Doxygen can't be found.")
  endif()
  #-- Configure the Template Doxyfile for our specific project
  configure_file(Doxyfile.in ${PROJECT_BINARY_DIR}/Doxyfile @ONLY IMMEDIATE)
  #-- Add a custom target to run Doxygen when ever the project is built
  add_custom_target (docs
  ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/Doxyfile SOURCES ${PROJECT_BINARY_DIR}/Doxyfile)
  # IF you do NOT want the documentation to be generated EVERY time you build the project
  # then leave out the 'ALL' keyword from the above command.
endif()

MESSAGE("The compiler is: ${CMAKE_CXX_COMPILER}")

################################################################################
# Subdirectories
add_subdirectory (./src/containers/sequential/tests)
add_subdirectory (./src/containers/hierarchical/tests)
add_subdirectory (./src/grid/)
add_subdirectory (./src/output/)
add_subdirectory (./src/solver/fv)
add_subdirectory (./src/geometry/)
################################################################################
