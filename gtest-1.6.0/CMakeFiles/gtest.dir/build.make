# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/Cellar/cmake/2.8.10.2/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/2.8.10.2/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/local/Cellar/cmake/2.8.10.2/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/gnzlbg/projects/sideprojects/HOM3/home

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/gnzlbg/projects/sideprojects/HOM3/home

# Include any dependencies generated for this target.
include gtest-1.6.0/CMakeFiles/gtest.dir/depend.make

# Include the progress variables for this target.
include gtest-1.6.0/CMakeFiles/gtest.dir/progress.make

# Include the compile flags for this target's objects.
include gtest-1.6.0/CMakeFiles/gtest.dir/flags.make

gtest-1.6.0/CMakeFiles/gtest.dir/src/gtest-all.cc.o: gtest-1.6.0/CMakeFiles/gtest.dir/flags.make
gtest-1.6.0/CMakeFiles/gtest.dir/src/gtest-all.cc.o: gtest-1.6.0/src/gtest-all.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/gnzlbg/projects/sideprojects/HOM3/home/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object gtest-1.6.0/CMakeFiles/gtest.dir/src/gtest-all.cc.o"
	cd /Users/gnzlbg/projects/sideprojects/HOM3/home/gtest-1.6.0 && /usr/local/bin/clang++   $(CXX_DEFINES) $(CXX_FLAGS) -arch x86_64 -Wall -Wextra -std=c++11 -stdlib=libc++ -I/Users/gnzlbg/src/libcxx/libcxx/include -pedantic -Wshadow -Woverloaded-virtual -pedantic-errors -Wcast-align -Wcomment -Wcast-qual  -Wchar-subscripts  -Wdisabled-optimization -Wfloat-equal -Wformat -Wformat=2 -Winvalid-pch -Wformat-nonliteral -Wformat-security -Wformat-y2k -Wimport  -Winit-self  -Winline -Wreturn-type  -Wmissing-braces -Wmissing-field-initializers -Wmissing-include-dirs  -Wredundant-decls -Wpacked -Wparentheses  -Wpointer-arith -Wsequence-point  -Wsign-compare  -Wstack-protector -Wstrict-aliasing -Wstrict-aliasing=2 -Wswitch -Wswitch-default  -Wtrigraphs -Wuninitialized -Wunknown-pragmas  -Wunreachable-code -Wunused -Wunused-function -Wunused-label  -Wunused-parameter -Wunused-value  -Wunused-variable -Wvariadic-macros -Wvolatile-register-var  -Wwrite-strings -Woverloaded-virtual -Wsign-promo -Wstrict-overflow=5 -Wswitch-default -DGTEST_USE_OWN_TR1_TUPLE=1 -fdiagnostics-show-template-tree -ftemplate-backtrace-limit=0 -Wno-attributes  -DGTEST_HAS_PTHREAD=1   -o CMakeFiles/gtest.dir/src/gtest-all.cc.o -c /Users/gnzlbg/projects/sideprojects/HOM3/home/gtest-1.6.0/src/gtest-all.cc

gtest-1.6.0/CMakeFiles/gtest.dir/src/gtest-all.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/gtest.dir/src/gtest-all.cc.i"
	cd /Users/gnzlbg/projects/sideprojects/HOM3/home/gtest-1.6.0 && /usr/local/bin/clang++  $(CXX_DEFINES) $(CXX_FLAGS) -arch x86_64 -Wall -Wextra -std=c++11 -stdlib=libc++ -I/Users/gnzlbg/src/libcxx/libcxx/include -pedantic -Wshadow -Woverloaded-virtual -pedantic-errors -Wcast-align -Wcomment -Wcast-qual  -Wchar-subscripts  -Wdisabled-optimization -Wfloat-equal -Wformat -Wformat=2 -Winvalid-pch -Wformat-nonliteral -Wformat-security -Wformat-y2k -Wimport  -Winit-self  -Winline -Wreturn-type  -Wmissing-braces -Wmissing-field-initializers -Wmissing-include-dirs  -Wredundant-decls -Wpacked -Wparentheses  -Wpointer-arith -Wsequence-point  -Wsign-compare  -Wstack-protector -Wstrict-aliasing -Wstrict-aliasing=2 -Wswitch -Wswitch-default  -Wtrigraphs -Wuninitialized -Wunknown-pragmas  -Wunreachable-code -Wunused -Wunused-function -Wunused-label  -Wunused-parameter -Wunused-value  -Wunused-variable -Wvariadic-macros -Wvolatile-register-var  -Wwrite-strings -Woverloaded-virtual -Wsign-promo -Wstrict-overflow=5 -Wswitch-default -DGTEST_USE_OWN_TR1_TUPLE=1 -fdiagnostics-show-template-tree -ftemplate-backtrace-limit=0 -Wno-attributes  -DGTEST_HAS_PTHREAD=1   -E /Users/gnzlbg/projects/sideprojects/HOM3/home/gtest-1.6.0/src/gtest-all.cc > CMakeFiles/gtest.dir/src/gtest-all.cc.i

gtest-1.6.0/CMakeFiles/gtest.dir/src/gtest-all.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/gtest.dir/src/gtest-all.cc.s"
	cd /Users/gnzlbg/projects/sideprojects/HOM3/home/gtest-1.6.0 && /usr/local/bin/clang++  $(CXX_DEFINES) $(CXX_FLAGS) -arch x86_64 -Wall -Wextra -std=c++11 -stdlib=libc++ -I/Users/gnzlbg/src/libcxx/libcxx/include -pedantic -Wshadow -Woverloaded-virtual -pedantic-errors -Wcast-align -Wcomment -Wcast-qual  -Wchar-subscripts  -Wdisabled-optimization -Wfloat-equal -Wformat -Wformat=2 -Winvalid-pch -Wformat-nonliteral -Wformat-security -Wformat-y2k -Wimport  -Winit-self  -Winline -Wreturn-type  -Wmissing-braces -Wmissing-field-initializers -Wmissing-include-dirs  -Wredundant-decls -Wpacked -Wparentheses  -Wpointer-arith -Wsequence-point  -Wsign-compare  -Wstack-protector -Wstrict-aliasing -Wstrict-aliasing=2 -Wswitch -Wswitch-default  -Wtrigraphs -Wuninitialized -Wunknown-pragmas  -Wunreachable-code -Wunused -Wunused-function -Wunused-label  -Wunused-parameter -Wunused-value  -Wunused-variable -Wvariadic-macros -Wvolatile-register-var  -Wwrite-strings -Woverloaded-virtual -Wsign-promo -Wstrict-overflow=5 -Wswitch-default -DGTEST_USE_OWN_TR1_TUPLE=1 -fdiagnostics-show-template-tree -ftemplate-backtrace-limit=0 -Wno-attributes  -DGTEST_HAS_PTHREAD=1   -S /Users/gnzlbg/projects/sideprojects/HOM3/home/gtest-1.6.0/src/gtest-all.cc -o CMakeFiles/gtest.dir/src/gtest-all.cc.s

gtest-1.6.0/CMakeFiles/gtest.dir/src/gtest-all.cc.o.requires:
.PHONY : gtest-1.6.0/CMakeFiles/gtest.dir/src/gtest-all.cc.o.requires

gtest-1.6.0/CMakeFiles/gtest.dir/src/gtest-all.cc.o.provides: gtest-1.6.0/CMakeFiles/gtest.dir/src/gtest-all.cc.o.requires
	$(MAKE) -f gtest-1.6.0/CMakeFiles/gtest.dir/build.make gtest-1.6.0/CMakeFiles/gtest.dir/src/gtest-all.cc.o.provides.build
.PHONY : gtest-1.6.0/CMakeFiles/gtest.dir/src/gtest-all.cc.o.provides

gtest-1.6.0/CMakeFiles/gtest.dir/src/gtest-all.cc.o.provides.build: gtest-1.6.0/CMakeFiles/gtest.dir/src/gtest-all.cc.o

# Object files for target gtest
gtest_OBJECTS = \
"CMakeFiles/gtest.dir/src/gtest-all.cc.o"

# External object files for target gtest
gtest_EXTERNAL_OBJECTS =

gtest-1.6.0/libgtest.a: gtest-1.6.0/CMakeFiles/gtest.dir/src/gtest-all.cc.o
gtest-1.6.0/libgtest.a: gtest-1.6.0/CMakeFiles/gtest.dir/build.make
gtest-1.6.0/libgtest.a: gtest-1.6.0/CMakeFiles/gtest.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX static library libgtest.a"
	cd /Users/gnzlbg/projects/sideprojects/HOM3/home/gtest-1.6.0 && $(CMAKE_COMMAND) -P CMakeFiles/gtest.dir/cmake_clean_target.cmake
	cd /Users/gnzlbg/projects/sideprojects/HOM3/home/gtest-1.6.0 && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/gtest.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
gtest-1.6.0/CMakeFiles/gtest.dir/build: gtest-1.6.0/libgtest.a
.PHONY : gtest-1.6.0/CMakeFiles/gtest.dir/build

gtest-1.6.0/CMakeFiles/gtest.dir/requires: gtest-1.6.0/CMakeFiles/gtest.dir/src/gtest-all.cc.o.requires
.PHONY : gtest-1.6.0/CMakeFiles/gtest.dir/requires

gtest-1.6.0/CMakeFiles/gtest.dir/clean:
	cd /Users/gnzlbg/projects/sideprojects/HOM3/home/gtest-1.6.0 && $(CMAKE_COMMAND) -P CMakeFiles/gtest.dir/cmake_clean.cmake
.PHONY : gtest-1.6.0/CMakeFiles/gtest.dir/clean

gtest-1.6.0/CMakeFiles/gtest.dir/depend:
	cd /Users/gnzlbg/projects/sideprojects/HOM3/home && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/gnzlbg/projects/sideprojects/HOM3/home /Users/gnzlbg/projects/sideprojects/HOM3/home/gtest-1.6.0 /Users/gnzlbg/projects/sideprojects/HOM3/home /Users/gnzlbg/projects/sideprojects/HOM3/home/gtest-1.6.0 /Users/gnzlbg/projects/sideprojects/HOM3/home/gtest-1.6.0/CMakeFiles/gtest.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : gtest-1.6.0/CMakeFiles/gtest.dir/depend

