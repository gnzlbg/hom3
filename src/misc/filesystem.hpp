#ifndef HOM3_MISC_FILESYSTEM_HPP_
#define HOM3_MISC_FILESYSTEM_HPP_
////////////////////////////////////////////////////////////////////////////////
/// \file \brief This file implements filesystem utilities.
////////////////////////////////////////////////////////////////////////////////
/// Includes:
#include <boost/filesystem.hpp>
#include "globals.hpp"
/// Options:
#define ENABLE_DBG_ 0
#include "misc/dbg.hpp"
////////////////////////////////////////////////////////////////////////////////
namespace hom3 { namespace io {
////////////////////////////////////////////////////////////////////////////////

/// \brief File related functionality
namespace file {

/// \brief Does file \p fileName exist?
bool exists(const String fileName, const boost::mpi::communicator& comm) {
  using namespace boost;
  return mpi::root_do(comm, [=]() { return filesystem::exists(fileName); });
}

/// \brief Removes the file \p fileName _if_ it exists.
/// \returns True if file existed, false otherwise.
bool remove(const String fileName, const boost::mpi::communicator& comm) {
  using namespace boost;
  return mpi::root_do(comm, [=]() { return filesystem::remove(fileName); });
}

}  // namespace file

////////////////////////////////////////////////////////////////////////////////
}  // namespace io
}  // namespace hom3
////////////////////////////////////////////////////////////////////////////////
#undef ENABLE_DBG_
////////////////////////////////////////////////////////////////////////////////
#endif
