#ifndef MATH_HPP_
#define MATH_HPP_
////////////////////////////////////////////////////////////////////////////////
/// \file \brief Mathematical functions
////////////////////////////////////////////////////////////////////////////////
/// Includes:
#include <limits>
#include "types.hpp"
////////////////////////////////////////////////////////////////////////////////
namespace hom3 {
////////////////////////////////////////////////////////////////////////////////

/// \brief Mathematical functions and helper classes
namespace math {

////////////////////////////////////////////////////////////////////////////////
// The following code is from the Google Test Testing Framework.
////////////////////////////////////////////////////////////////////////////////
// Copyright 2005, Google Inc.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
// copyright notice, this list of conditions and the following disclaimer
// in the documentation and/or other materials provided with the
// distribution.
//     * Neither the name of Google Inc. nor the names of its
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Authors: wan@google.com (Zhanyong Wan), eefacm@gmail.com (Sean Mcafee)
//
// The Google C++ Testing Framework (Google Test)


/// \brief This template class serves as a compile-time function from size to
/// type.  It maps a size in bytes to a primitive type with that
/// size. e.g.
///
/// \code
///   TypeWithSize<4>::UInt
/// \endcode
///
/// is typedef-ed to be unsigned int (unsigned integer made up of 4
/// bytes).
///
/// Such functionality should belong to STL, but I cannot find it
/// there.
///
/// Google Test uses this class in the implementation of floating-point
/// comparison.
///
/// For now it only handles UInt (unsigned int) as that's all Google Test
/// needs.  Other types can be easily added in the future if need
/// arises.
template <size_t size> struct TypeWithSize {
  /// This prevents the user from using TypeWithSize<N> with incorrect
  /// values of N.
  typedef void UInt;
};

/// The specialization for size 4.
template <> struct TypeWithSize<4> {
  /// unsigned int has size 4 in both gcc and MSVC.
  ///
  /// As base/basictypes.h doesn't compile on Windows, we cannot use
  /// uint32, uint64, and etc here.
  typedef int Int;
  typedef unsigned int UInt;
};

/// The specialization for size 8.
template <> struct TypeWithSize<8> {
  typedef long long Int;  // NOLINT
  typedef unsigned long long UInt;  // NOLINT
};


/// \brief This template class represents an IEEE floating-point number
/// (either single-precision or double-precision, depending on the
/// template parameters).
///
/// The purpose of this class is to do more sophisticated number
/// comparison.  (Due to round-off error, etc, it's very unlikely that
/// two floating-points will be equal exactly.  Hence a naive
/// comparison by the == operation often doesn't work.)
///
/// Format of IEEE floating-point:
///
///   The most-significant bit being the leftmost, an IEEE
///   floating-point looks like
///
///     sign_bit exponent_bits fraction_bits
///
///   Here, sign_bit is a single bit that designates the sign of the
///   number.
///
///   For float, there are 8 exponent bits and 23 fraction bits.
///
///   For double, there are 11 exponent bits and 52 fraction bits.
///
///   More details can be found at
///   http://en.wikipedia.org/wiki/IEEE_floating-point_standard.
///
/// Template parameter:
///
/// \param [T] RawType: the raw floating-point type (either float or double)
template <typename RawType> struct FloatingPoint {
  /// Defines the unsigned integer type that has the same size as the
  /// floating point number.
  using Bits = typename TypeWithSize<sizeof(RawType)>::UInt;

  /// Constants.

  /// # of bits in a number.
  static const constexpr size_t kBitCount = 8*sizeof(RawType);

  /// # of fraction bits in a number.
  static const constexpr size_t kFractionBitCount
  = std::numeric_limits<RawType>::digits - 1;

  /// # of exponent bits in a number.
  static const constexpr size_t kExponentBitCount
  = kBitCount - 1 - kFractionBitCount;

  /// The mask for the sign bit.
  static const constexpr Bits kSignBitMask
  = static_cast<Bits>(1) << (kBitCount - 1);

  // The mask for the fraction bits.
  static const constexpr Bits kFractionBitMask
  = (~static_cast<Bits>(0) >> (kExponentBitCount + 1));

  /// The mask for the exponent bits.
  static const constexpr Bits kExponentBitMask
  = ~(kSignBitMask | kFractionBitMask);

  /// How many ULP's (Units in the Last Place) we want to tolerate when
  /// comparing two numbers.  The larger the value, the more error we
  /// allow.  A 0 value means that two numbers must be exactly the same
  /// to be considered equal.
  ///
  /// The maximum error of a single floating-point operation is 0.5
  /// units in the last place.  On Intel CPU's, all floating-point
  /// calculations are done with 80-bit precision, while double has 64
  /// bits.  Therefore, 4 should be enough for ordinary use.
  ///
  /// See the following article for more details on ULP:
  /// http://www.cygnus-software.com/papers/comparingfloats/comparingfloats.htm.
  static const constexpr size_t kMaxUlps = 4;

  /// Constructs a FloatingPoint from a raw floating-point number.
  ///
  /// On an Intel CPU, passing a non-normalized NAN (Not a Number)
  /// around may change its bits, although the new value is guaranteed
  /// to be also a NAN.  Therefore, don't expect this constructor to
  /// preserve the bits in x when x is a NAN.
  explicit FloatingPoint(const RawType& x) noexcept { u_.value_ = x; }

  /// Static methods

  /// Reinterprets a bit pattern as a floating-point number.
  ///
  /// This function is needed to test the AlmostEquals() method.
  static inline constexpr RawType ReinterpretBits(const Bits bits) noexcept {
    FloatingPoint fp(0);
    fp.u_.bits_ = bits;
    return fp.u_.value_;
  }

  /// Returns the floating-point number that represent positive infinity.
  static inline constexpr RawType Infinity() noexcept {
    return ReinterpretBits(kExponentBitMask);
  }

  /// Non-static methods

  /// Returns the bits that represents this number.
  inline constexpr const Bits &bits() const noexcept { return u_.bits_; }

  /// Returns the exponent bits of this number.
  inline constexpr Bits exponent_bits() const noexcept
  { return kExponentBitMask & u_.bits_; }

  /// Returns the fraction bits of this number.
  inline constexpr Bits fraction_bits() const noexcept
  { return kFractionBitMask & u_.bits_; }

  /// Returns the sign bit of this number.
  inline constexpr Bits sign_bit() const noexcept
  { return kSignBitMask & u_.bits_; }

  /// Returns true iff this is NAN (not a number).
  inline constexpr bool is_nan() const noexcept {
    /// It's a NAN if the exponent bits are all ones and the fraction
    /// bits are not entirely zeros.
    return (exponent_bits() == kExponentBitMask) && (fraction_bits() != 0);
  }

  /// Returns true iff this number is at most kMaxUlps ULP's away from
  /// rhs.  In particular, this function:
  ///
  ///   - returns false if either number is (or both are) NAN.
  ///   - treats really large numbers as almost equal to infinity.
  ///   - thinks +0.0 and -0.0 are 0 DLP's apart.
  inline constexpr bool AlmostEquals(const FloatingPoint& rhs) const noexcept {
    /// The IEEE standard says that any comparison operation involving
    /// a NAN must return false.
    if (is_nan() || rhs.is_nan()) return false;

    return DistanceBetweenSignAndMagnitudeNumbers(u_.bits_, rhs.u_.bits_)
        <= kMaxUlps;
  }

 private:
  /// The data type used to store the actual floating-point number.
  union FloatingPointUnion {
    RawType value_;  ///< The raw floating-point number.
    Bits bits_;      ///< The bits that represent the number.
  };

  /// \brief Converts an integer from the sign-and-magnitude representation to
  /// the biased representation.  More precisely, let N be 2 to the
  /// power of (kBitCount - 1), an integer x is represented by the
  /// unsigned number x + N.
  ///
  /// For instance,
  ///
  ///   -N + 1 (the most negative number representable using
  ///          sign-and-magnitude) is represented by 1;
  ///   0      is represented by N; and
  ///   N - 1  (the biggest number representable using
  ///          sign-and-magnitude) is represented by 2N - 1.
  ///
  /// Read http://en.wikipedia.org/wiki/Signed_number_representations
  /// for more details on signed number representations.
  inline constexpr static Bits SignAndMagnitudeToBiased
  (const Bits &sam) noexcept {
    if (kSignBitMask & sam) {
      /// sam represents a negative number.
      return ~sam + 1;
    } else {
      /// sam represents a positive number.
      return kSignBitMask | sam;
    }
  }

  /// \brief Given two numbers in the sign-and-magnitude representation,
  /// returns the distance between them as an unsigned number.
  inline constexpr static Bits DistanceBetweenSignAndMagnitudeNumbers
  (const Bits &sam1, const Bits &sam2) noexcept {
    const Bits biased1 = SignAndMagnitudeToBiased(sam1);
    const Bits biased2 = SignAndMagnitudeToBiased(sam2);
    return (biased1 >= biased2) ? (biased1 - biased2) : (biased2 - biased1);
  }

  FloatingPointUnion u_;
};
////////////////////////////////////////////////////////////////////////////////

/// \brief Creates a floating point type
///
/// Note: this is an alias for the GoogleTest FPA class
///
/// Example:
/// \code
/// if ( FP(a) == FP(b) ) {...}
/// \endcode
///
template<class T> using FP = FloatingPoint<T>;

/// \todo Compares two floating point values for equality
template<class T>
static inline constexpr bool approx(const T& a, const T& b) noexcept
{ return FP<T>(a).AlmostEquals(FP<T>(b)); }

/// Takes the difference of two unsigned integers (without wrapping)
template<class T>
static inline constexpr T absdiff(const T& a, const T& b) noexcept
{ return (a > b) ? (a - b) : (b - a); }

template <class T>
static inline constexpr int signum(const T& x, std::false_type) noexcept
{ return T{0} < x; }

template <class T>
static inline constexpr int signum(const T& x, std::true_type) noexcept
{ return (T{0} < x) - (x < T{0}); }

/// \brief Returns the signum of \p x, i.e. -1 if x < 0, +1 if x > 0, and 0
/// otherwise
template <class T>
static inline constexpr int signum(const T& x) noexcept
{ return signum(x, std::is_signed<T>()); }

/// Useful traits
namespace traits {

template<class T> struct is_long_int : std::false_type {};
template<> struct is_long_int<Ind> : std::true_type {};
template<> struct is_long_int<Int> : std::true_type {};

} // namespace traits

/// \brief Compile-time numeric functions
namespace ct {

/// \brief Computes the logarithm log_b(n) for (b,n) integers
///
/// Note: Overflows for N > 2^64.
constexpr long ilog(const long b, const long n) {
  return n == 0l ? 0l :
         n == 1l ? 0l :
        1l + ilog(b,n/b);
}

/// \brief Computes b^e for (b,e) integers
///
/// Note: overflows for e >= (64-1) / ilog(b,2);
constexpr long ipow(const long b, const long e){
  return e == 0l ? 1l :
         b * ipow(b,e-1l);
}

} // namespace ct

/// \brief Run-time numeric functions
namespace rt {

namespace detail_ {

using math::traits::is_long_int;

template<class Integral1, class Integral2>
Integral1 ilog(Integral1 b, Integral2 n){
  return n==0 ? 0 : ( n==1 ? 0 : 1 + ilog(b,n/b) ); // used only in asserts
}

template<class Integral1, class Integral2>
bool ipow_overflow_long(const Integral1 b, const Integral2 e) {
  return e >= (64-1) / ilog(b, static_cast<Integral1>(2));
}

template<class Integral1, class Integral2>
bool ipow_overflow_not_long(const Integral1 b, const Integral2 e) {
  return e >= (32-1) / ilog(b, static_cast<Integral1>(2));
}

template<class Integral1, class Integral2>
bool ipow_overflow(Integral1 b, Integral2 e) {
  return is_long_int<Integral1>::value
      ? ipow_overflow_long(b,e)
      : ipow_overflow_not_long(b,e);
}

} // detail_ namespace

/// Computes integer pow using exponentiation by squaring
/// Complexiy O(log(e))
template<class Integral1, class Integral2>
Integral1 ipow(Integral1 b, Integral2 e) {
  ASSERT(!detail_::ipow_overflow(b,e),"This expression overflows");
  Integral1 result = 1;
  while (e) {
    if (e & 1) { result *= b; }
    e >>= 1;
    b *= b;
  }
  return result;
}

} // namespace rt

} // namespace math

////////////////////////////////////////////////////////////////////////////////
} // namespace hom3
////////////////////////////////////////////////////////////////////////////////
#endif
