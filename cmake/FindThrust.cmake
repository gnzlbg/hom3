# - Try to find thrust
# Once done this will define
# Thrust_FOUND - System has LibXml2

find_package(PkgConfig)

find_path(THRUST_INCLUDE_DIR
    HINTS /usr/include/cuda /usr/local/include
    NAMES thrust/version.h
    PATHS /usr/include/cuda /usr/local/include
    PATH_SUFFIXES thrust
    DOC "Thrust headers"
)

include(FindPackageHandleStandardArgs)
mark_as_advanced(THRUST_INCLUDE_DIR)
