# - Try to find hdf5 lite
# Once done this will define
#  HDF5LITE_FOUND - System has LibXml2

find_package(PkgConfig)

find_path(HDF5LITE_INCLUDE_DIR hdf5_hl.h
          HINTS ${PC_HDF5LITE_INCLUDEDIR})

find_library(HDF5LITE_LIBRARY NAMES hdf5_hl
             HINTS ${PC_HDF5LITE_LIBDIR})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(HDF5LITE DEFAULT_MSG
                                  HDF5LITE_LIBRARY HDF5LITE_INCLUDE_DIR)

mark_as_advanced(HDF5LITE_INCLUDE_DIR HDF5LITE_LIBRARY)
