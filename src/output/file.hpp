#ifndef HOM3_FILE_HPP_
#define HOM3_FILE_HPP_
////////////////////////////////////////////////////////////////////////////////
/// \file \brief This file includes the selected I/O system
///
/// \file \long This file includes the I/O system. To select:
///   - PNetCDF I/O: define the HOM3_ENABLE_IO_PNETCDF macro
///   - HDF5 I/O: define the HOM3_ENABLE_IO_HDF5 macro
///
////////////////////////////////////////////////////////////////////////////////
/// Includes:
#include "../globals.hpp"
#include "modes.hpp"
#include "data_types.hpp"
#include "information.hpp"
#include "filesystem.hpp"
////////////////////////////////////////////////////////////////////////////////

#ifdef HOM3_ENABLE_IO_PNETCDF
  #include "pnetcdf.hpp"
#elseif HOM3_ENABLE_IO_HDF5
  #include "hdf5.hpp"
#else
  #pragma error "unknown IO system"
#endif

////////////////////////////////////////////////////////////////////////////////
#endif
