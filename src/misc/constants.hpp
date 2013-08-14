#ifndef HOM3_MISC_CONSTANTS_HPP_
#define HOM3_MISC_CONSTANTS_HPP_
////////////////////////////////////////////////////////////////////////////////
/// \file \brief Includes numeric constants
////////////////////////////////////////////////////////////////////////////////
/// Includes:
#include <limits>
#include <boost/math/constants/constants.hpp>
#include "types.hpp"
#include "traits.hpp"
////////////////////////////////////////////////////////////////////////////////

namespace hom3 {

/// Sometimes you want to represent an empty value,
/// examples:
/// - no nghbr has been found
/// - node id does not exist
///
/// However, node ids might have type unsigned long long, and neighbor positions
/// might have type unsigned int. That is you cannot fix this unknown value to
/// either max<ulonglong> or max<uint> because one might be too large for the
/// other. The following functions try to fix this. This solution is better
/// than what i previously had. For integers, using a Maybe type like
/// std/boost::optional seems overkill.
///
/// \todo constrain these template s.t. they only work for hom3 types
///
/// They are in numeric constants because they define a family of numeric
/// constants and operations on them, and because they were previously constants
/// and were already here.

template<class T, class B> struct Integer;

namespace constants { namespace detail {

/// \brief \returns invalid value for identifier of type T
template<class T> constexpr T invalid_value_(T) {
  return std::numeric_limits<T>::max();
}

constexpr Num invalid_value_(Num) {
  static_assert(std::numeric_limits<Num>::has_quiet_NaN, "no NaN support!");
  return std::numeric_limits<Num>::quiet_NaN();
}

template<class T,class B> constexpr Integer<T,B> invalid_value_(Integer<T,B>) {
  return Integer<T,B>{invalid_value_(T())};
}

}} // constants::detail namespace

template<class T> constexpr T invalid() { return constants::detail::invalid_value_(T()); }

/// \brief \returns [bool] Has the identifier "o" a valid value?
template<class T> static constexpr bool is_valid(const T& o) {
  return o != invalid<T>();
}

namespace math {
static constexpr Num eps = std::numeric_limits<Num>::epsilon();
static constexpr Num pi = boost::math::constants::pi<Num>();

template<class T, EnableIf<traits::is_eigen_matrix<T>> = traits::dummy>
constexpr auto zero(T&&) -> decltype(T::Zero()) {
  return T::Zero();
}
template<class T, DisableIf<traits::is_eigen_matrix<T>> = traits::dummy>
constexpr auto zero(T&&) -> decltype(T(0)) {
  return T(0);
}

template<class T> constexpr auto dimensions(T&& t) -> decltype(t.cols()) {
  return t.cols();
}

} // math

////////////////////////////////////////////////////////////////////////////////
} // hom3 namespace
////////////////////////////////////////////////////////////////////////////////
#endif
