# - Try to find pnetcdf
# Once done this will define
#  PNETCDF_FOUND - System has LibXml2

find_package(PkgConfig)

find_path(PNETCDF_INCLUDE_DIR pnetcdf.h
          HINTS ${PC_PNETCDF_INCLUDEDIR})

find_library(PNETCDF_LIBRARY NAMES pnetcdf
             HINTS ${PC_PNETCDF_LIBDIR})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PNETCDF DEFAULT_MSG
                                  PNETCDF_LIBRARY PNETCDF_INCLUDE_DIR)

mark_as_advanced(PNETCDF_INCLUDE_DIR PNETCDF_LIBRARY)
