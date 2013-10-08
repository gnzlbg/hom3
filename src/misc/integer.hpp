#ifndef HOM3_MISC_INTEGER_HPP_
#define HOM3_MISC_INTEGER_HPP_
////////////////////////////////////////////////////////////////////////////////
#include <boost/iterator/counting_iterator.hpp>
#include <type_traits>
#include <limits>
#include <string>
#include <algorithm>
#include "constants.hpp"
////////////////////////////////////////////////////////////////////////////////
namespace hom3 {
////////////////////////////////////////////////////////////////////////////////

/// \name Index types
///@{
/// \brief Implements an integer type
template<class T, class B = void> struct Integer {
  using value_type = T;

  /// \name Asignment operators
  ///@{
  constexpr Integer() noexcept : value(invalid<T>()) {}
  constexpr Integer(const Integer& other) noexcept : value(other.value) {}
  constexpr Integer(Integer&& other) noexcept : value(other.value) {}
  constexpr explicit Integer(const T& other) noexcept : value(other) {}

  constexpr inline Integer& operator=(const Integer& other) noexcept {
    value = other.value;
    return *this;
  }
  constexpr inline Integer& operator=(Integer&& other) noexcept {
    value = other.value;
    return *this;
  }
  constexpr inline Integer& operator=(const T& other) noexcept {
    value = other;
    return *this;
  }
  ///@}

  /// \name Conversion operators
  ///@{

  explicit constexpr inline operator       T()       noexcept { return value; }
  explicit constexpr inline operator const T() const noexcept { return value; }

  ///@}

  /// \name Compound assignment +=, -=, *=, /=
  ///@{
  constexpr inline Integer& operator+=(const Integer& other) noexcept {
    value += other.value;
    return *this;
  }
  constexpr inline Integer& operator-=(const Integer& other) noexcept {
    value -= other.value;
    return *this;
  }
  constexpr inline Integer& operator*=(const Integer& other) noexcept {
    value *= other.value;
    return *this;
  }
  constexpr inline Integer& operator/=(const Integer& other) noexcept {
    value /= other.value;
    return *this;
  }
  ///@}

  /// \name Arithmetic operators +,-,*,/,unary -
  ///@{
  constexpr friend inline
  Integer operator+(Integer a, const Integer& b) noexcept { return a += b; }
  constexpr friend inline
  Integer operator-(Integer a, const Integer& b) noexcept { return a -= b; }
  constexpr friend inline
  Integer operator*(Integer a, const Integer& b) noexcept { return a *= b; }
  constexpr friend inline
  Integer operator/(Integer a, const Integer& b) noexcept { return a /= b; }

  constexpr inline Integer operator-() noexcept {
    static_assert(std::is_signed<T>::value, "Can't negate an unsigned type!");
    return Integer{-value};
  }
  ///@}

  /// \name Prefix increment operators ++(),--()
  ///@{
  constexpr inline Integer& operator++() noexcept {
    ++value;
    return *this;
  }
  constexpr inline Integer& operator--() noexcept {
    --value;
    return *this;
  }
  ///@}

  /// \name Postfix increment operators ()++,()--
  ///@{
  constexpr inline Integer operator++(int) noexcept {
    Integer tmp(*this);
    ++(*this);
    return tmp;
  }
  constexpr inline Integer operator--(int) noexcept {
    Integer tmp(*this);
    --(*this);
    return tmp;
  }
  ///@}

  /// \name Comparison operators ==, !=, <, >, <=, >=
  ///@{
  constexpr friend inline
  bool operator==(const Integer& a, const Integer& b) noexcept
  { return a.value == b.value; }
  constexpr friend inline
  bool operator<=(const Integer& a, const Integer& b) noexcept
  { return a.value <= b.value; }
  constexpr friend inline
  bool operator<(const Integer& a, const Integer& b)  noexcept
  { return a.value < b.value; }  // return a <= b && !(a == b) -> slower?
  constexpr friend inline
  bool operator!=(const Integer& a, const Integer& b) noexcept
  { return !(a == b); }
  constexpr friend inline
  bool operator>(const Integer& a, const Integer& b)  noexcept
  { return !(a <= b); }
  constexpr friend inline
  bool operator>=(const Integer& a, const Integer& b) noexcept
  { return !(a < b); }
  ///@}

  /// \brief swap
  constexpr friend inline void swap(Integer&& a, Integer&& b) noexcept {
    using std::swap;
    swap(a.value, b.value);
  }

  /// \brief to_string
  friend inline std::string to_string(const Integer a)  {
    return std::to_string(a.value);
  }

  /// \name Access operator
  ///@{
  constexpr inline T& operator()()       noexcept { return value; }
  constexpr inline T  operator()() const noexcept { return value; }
  ///@}

  T value;

  template<class C, class CT>
  friend inline std::basic_ostream<C, CT>& operator<<
  (std::basic_ostream<C, CT>& o, const Integer<T, B>& i) {
    return o << i();
  }
};
///@}

template<class T> struct is_integer_
{ static const bool value = false; };
template<class T, class U> struct is_integer_<Integer<T, U>>
{ static const bool value = true; };


template<class T> struct is_integer
{ static const bool value = is_integer_<
  std::remove_reference_t<std::remove_cv_t<T>>>::value;
};


template<class T, EnableIf<is_integer<T>> = traits::dummy>
constexpr inline auto&& primitive_cast(T&& t) {
  static_assert(is_integer<T>::value, "T must be an integer!");
  return t();
}

template<class T, EnableIf<is_integer<T>> = traits::dummy>
constexpr inline auto primitive_cast(const T& t) {
  static_assert(is_integer<T>::value, "T must be an integer!");
  return t();
}

template<class T, EnableIf<is_integer<T>> = traits::dummy>
constexpr inline auto& primitive_cast(T& t) {
  static_assert(is_integer<T>::value, "T must be an integer!");
  return t();
}


template<class T, DisableIf<is_integer<T>> = traits::dummy>
constexpr inline auto&& primitive_cast(T&& t) {
  static_assert(!is_integer<T>::value, "T can't be an integer!");
  return std::forward<T>(t);
}

template<class T, DisableIf<is_integer<T>> = traits::dummy>
constexpr inline auto primitive_cast(const T& t) {
  static_assert(!is_integer<T>::value, "T can't be an integer!");
  return t;
}

template<class T, DisableIf<is_integer<T>> = traits::dummy>
constexpr inline auto& primitive_cast(T& t) {
  static_assert(!is_integer<T>::value, "T can't be an integer!");
  return t;
}


////////////////////////////////////////////////////////////////////////////////
}  // namespace hom3
////////////////////////////////////////////////////////////////////////////////

// this feels wrong:
namespace std {
template<class T, class B>
class numeric_limits<hom3::Integer<T, B>> : public numeric_limits<T> {
 public:
  static const bool is_specialized = true;
};
}  // namespace std

////////////////////////////////////////////////////////////////////////////////
#endif
