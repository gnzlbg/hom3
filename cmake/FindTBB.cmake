# Find Intel's TBB library: https://www.threadingbuildingblocks.org/

find_package(PkgConfig)

find_path(TBB_INCLUDE_DIR tbb/tbb.h
  HINTS /usr/local/include ${PC_TBB_INCLUDE_DIR}
)

find_library(TBB_CORE_LIBRARY NAMES tbb)
find_library(TBB_MALLOC_LIBRARY NAMES tbbmalloc)


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(TBB  DEFAULT_MSG  TBB_CORE_LIBRARY  TBB_INCLUDE_DIR)

mark_as_advanced(TBB_INCLUDE_DIR TBB_LIBRARY)
