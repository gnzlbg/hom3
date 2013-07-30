#ifndef HOM3_HDF5_FILE_
#define HOM3_HDF5_FILE_
////////////////////////////////////////////////////////////////////////////////
/// Include:
#include <hdf5.h>
#include <boost/optional.hpp>

////////////////////////////////////////////////////////////////////////////////
namespace hdf5 {
////////////////////////////////////////////////////////////////////////////////

struct File {

  File(const std::string fName) : fileName(fName) {
    open_excl(fileName);
  }
  ~File() {
    close(file_id);
  }

  void open_excl(const std::string fName) {
    // Create a new file using default properties:
    file_id = H5Fcreate(fName, H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
  }

  void open(const std::string fName) {
    // Create a new file using default properties:
    file_id = H5Fcreate(fName, H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
  }

  void close() {
    if(file_id) {
      status = H5Fclose(file_id);
    }
  }



 private:
  boost::optional<hid_t> file_id;
  herr_t status;
};






////////////////////////////////////////////////////////////////////////////////
} // namespace hdf5
////////////////////////////////////////////////////////////////////////////////
#endif
