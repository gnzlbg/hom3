/// Includes:
#include "../globals.hpp"
#include "file.hpp"
/// External Includes:
#include "gtest/gtest.h"
/// Options:
#define ENABLE_DBG_ 0
#include "../misc/dbg.hpp"
////////////////////////////////////////////////////////////////////////////////

int argc = 0;
char** argv = NULL;
boost::mpi::environment env(argc,argv);


/// \test Test file handler constructors
TEST(file_test, constructor) {

  EXPECT_DEATH(
      io::File<io::file::type::PNetCDF>("does_not_exist.Netcdf", io::file::op::read()),
      "[\\S\\s]+");

}

/// \test Write/read data to file
TEST(file_test, write_read_data_to_file) {

  /// Write data to a file
  {
    io::File<io::file::type::PNetCDF>("test.Netcdf", io::file::op::write());

  }

  /// Read data from the file
  {


  }
}



////////////////////////////////////////////////////////////////////////////////
#undef ENABLE_DBG_
////////////////////////////////////////////////////////////////////////////////
