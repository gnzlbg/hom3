#ifndef HOM3_MISC_ERROR_HPP_
#define HOM3_MISC_ERROR_HPP_
////////////////////////////////////////////////////////////////////////////////
/// \file \brief Error-Handling functions
////////////////////////////////////////////////////////////////////////////////
/// INCLUDES:
#include <iostream>
#include <string>
#include <stdexcept>
#include "misc/functions.hpp"
////////////////////////////////////////////////////////////////////////////////
namespace hom3 {
////////////////////////////////////////////////////////////////////////////////

/// \brief Contains error-handling functions
namespace error {

/// \brief Throws a runtime error exception with the description given
/// in \p error_string.
std::runtime_error exception(std::string error_string) {
  std::runtime_error tmp("EXCEPTION: " + error_string);
  return tmp;
}

/// \brief Terminates the program in case of error
///
/// \warning this function never returns
[[noreturn, noinline, cold]] void terminate
(const std::string error_string, const std::string at) {
  std::cerr << "FATAL ERROR \"" << error_string << "\" (" << at << ")."
            << std::endl;
  exit(1);
}

}  // namespace error

/// \brief Terminates the program
#define TERMINATE(error_string)                 \
  hom3::error::terminate(error_string, AT_)     \

////////////////////////////////////////////////////////////////////////////////
}  // namespace hom3
////////////////////////////////////////////////////////////////////////////////
#endif
