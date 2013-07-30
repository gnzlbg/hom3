#ifndef HOM3_MISC_TRATIS_HPP_
#define HOM3_MISC_TRATIS_HPP_
////////////////////////////////////////////////////////////////////////////////
namespace traits {

namespace detail { enum class enabler {}; }
constexpr detail::enabler dummy = {};

template<class C, class T = detail::enabler>
using EnableIf = typename std::enable_if<C::value, T>::type;

template<class C, class T = detail::enabler>
using DisableIf = typename std::enable_if<!C::value, T>::type;

template<class T, T a, T b, class C = T> struct equal : std::integral_constant<bool, a == b> {};

struct lazy {};
struct strict {};


/// \brief Detects if T is an Eigen::Matrix type, in general it is not:
template<class T> struct is_eigen_matrix { static const bool value = false; };

///\brief Detects if T is an Eigen::Matrix type, sometimes it is:
template<class _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
struct is_eigen_matrix<Eigen::Matrix<_Scalar,_Rows,_Cols,_Options,_MaxRows,_MaxCols>> {
  static const bool value = true;
};

} // traits namespace

using traits::lazy;
using traits::strict;
////////////////////////////////////////////////////////////////////////////////
#endif
