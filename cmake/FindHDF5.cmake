# - Try to find pnetcdf
# Once done this will define
#  HDF5_FOUND - System has LibXml2

find_package(PkgConfig)

find_path(HDF5_INCLUDE_DIR hdf5.h
          HINTS ${PC_HDF5_INCLUDEDIR})

find_library(HDF5_LIBRARY NAMES hdf5
             HINTS ${PC_HDF5_LIBDIR})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(HDF5 DEFAULT_MSG
                                  HDF5_LIBRARY HDF5_INCLUDE_DIR)

mark_as_advanced(HDF5_INCLUDE_DIR HDF5_LIBRARY)
