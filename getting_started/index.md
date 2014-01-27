---
layout: default
---

# Getting started 

Hom3 requires a fairly modern version of [clang](http://clang.llvm.org/)
(git-HEAD) and [libc++](http://libcxx.llvm.org/) as well as the latest version
of [cmake](http://www.cmake.org/), [boost](http://www.boost.org/) (in
particular: Range, MPI, Filesystem, Units, and their dependencies),
[eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page),
[hdf5](http://www.hdfgroup.org/HDF5/),
[doxygen](http://www.stack.nl/~dimitri/doxygen/),
[tbb](https://www.threadingbuildingblocks.org/),
[thrust](https://github.com/thrust/thrust), and
[mpich](http://www.mpich.org/). The latest
[gtest](https://code.google.com/p/googletest/) version is automatically fetched
from svn (an internet connection is required).

    ./configure.sh -d && make && ctest # should configure, compile, and run all hom3 tests

If configure fails you might need to tweak the `configure.sh` script to tell
cmake about your local library paths.

The configure script can set up debug (`-d`) and release (`-r`) builds as well
as sanitizer builds: `-a` for [address
sanitizer](http://clang.llvm.org/docs/AddressSanitizer.html), `-m` for [memory
sanitizer](http://clang.llvm.org/docs/MemorySanitizer.html), and `-t` for
[thread sanitizer](http://clang.llvm.org/docs/ThreadSanitizer.html). For more
help, see: `./configure.sh -h`.

To generate the documentation: `make docs`.
To check for style-guide conformance: `make style`.
