#ifndef HOM3_MISC_ASSERT_HPP_
#define HOM3_MISC_ASSERT_HPP_
////////////////////////////////////////////////////////////////////////////////
/// Includes:
#include<iostream>
////////////////////////////////////////////////////////////////////////////////

/// \brief Returns a string containing the position where it is defined. That
/// is, " at function _function_name_ in file _file_name_ at line _line_number_
#define AT_ (std::string(" at function ") + __PRETTY_FUNCTION__ + " in file " \
             + __FILE__ + " at line " + std::to_string(__LINE__))

/// \brief Asserts a condition in DEBUG mode. Useful for checking
/// pre/postconditions and invariants.
///
/// \warning if NDEBUG is defined the condition in the if clause is always
/// false and ASSERT compiles to nothing! Still this happens after the
/// compiler has verified the correctness of the ASSERT code.
#ifndef NDEBUG
#   define ASSERT(condition, message)           \
  do {                                          \
    if (!(condition)) {                                                 \
      std::cerr << "Assertion `" #condition "` failed in " << AT_       \
                << ": " << message << std::endl;                        \
          assert(false);                                                \
          std::exit(EXIT_FAILURE);                                      \
    }                                                                   \
  } while (false)
#else
#   define ASSERT(condition, message)                                   \
  do {                                                                  \
    if (0 && (condition)) {                                             \
      std::cerr << "Assertion `" #condition "` failed in " << AT_       \
                << ": " << message << std::endl;                        \
          std::exit(EXIT_FAILURE);                                      \
    }                                                                   \
  } while (false)
#endif

////////////////////////////////////////////////////////////////////////////////
#endif
