#ifndef HOM3_MISC_FUNCTIONS_HPP_
#define HOM3_MISC_FUNCTIONS_HPP_
////////////////////////////////////////////////////////////////////////////////
/// \file \brief Contains helper functions
////////////////////////////////////////////////////////////////////////////////
/// Includes:
#include <memory>
////////////////////////////////////////////////////////////////////////////////

/// \brief Performs a cold computation preventing inlining and branch prediction
template<class F>
[[noinline,cold]] auto cold_do(F&& f) -> decltype(f()) {
  return f();
}

/// \brief Creates a std::unique_ptr of type T (analog to make_shared)
/// \todo move into some memory.hpp header
template<class T, class... Args>
std::unique_ptr<T> make_unique(Args&&... args) {
  return std::unique_ptr<T>( new T(std::forward<Args>(args)...) );
}

////////////////////////////////////////////////////////////////////////////////
#endif
