#ifndef HOM3_CONTAINER_SEQUENTIAL_TEST_FIXEDCELL_HPP_
#define HOM3_CONTAINER_SEQUENTIAL_TEST_FIXEDCELL_HPP_
////////////////////////////////////////////////////////////////////////////////
/// \file \brief Fixed container for testing pourposes.
////////////////////////////////////////////////////////////////////////////////
// Includes:
#include <algorithm>
#include "containers/sequential.hpp"
/// Options:
#define ENABLE_DBG_ 0
#include "misc/dbg.hpp"
////////////////////////////////////////////////////////////////////////////////
namespace hom3 { namespace container { namespace sequential {
////////////////////////////////////////////////////////////////////////////////

namespace test {

/// Fixed-container forward declarations:
namespace fixed {
template<SInd nd> struct Container;
template<SInd nd> struct Reference;
template<SInd nd> struct Value;
}  // namespace fixed

}  // namspace test

/// Fixed Test Container trait:
template<SInd nd> struct traits<test::fixed::Container<nd>> {
  using container_type  = tag::fixed_nodes;
  using value_type      = typename test::fixed::Value<nd>;
  using reference       = typename test::fixed::Reference<nd>;
  using cell_index_type = Ind;
};

namespace test { namespace fixed {

/// \brief Fixed-size sequential test container
template<SInd nd_> struct Container : container::Sequential<Container<nd_>> {
  static const SInd nd = nd_;
  friend container::Sequential<Container<nd>>;

  /// Aliases:
  template<template <SInd> class V_, SInd nd__ = 1> using M = container::Matrix
  <Container, container::matrix::tag::Cell, V_, Ind, SInd, nd__>;

  Container(const Ind n) noexcept
    : container::Sequential<Container<nd>>(n, "fixed_container")
    ,  mInt(this, "mInt" )
    ,  mNum(this, "mNum" )
    , mIntA(this, "mIntA")
    , mNumA(this, "mNumA")
  {}

  /// Data:
  M<IntM>     mInt;
  M<NumM>     mNum;
  M<IntM, nd> mIntA;
  M<NumM, nd> mNumA;

 private:
  /// \brief Resets variables of cell \p cIdx
  inline void reset_cell(const Ind cIdx) noexcept {
    mInt(cIdx) = 0;
    mNum(cIdx) = 0;
    for (SInd d = 0; d < nd; ++d) {
      mIntA(cIdx, d) = 0;
      mNumA(cIdx, d) = 0;
    }
  }

  /// \brief Copies cell variables from cell \p fromIdx to \p toIdx
  inline void copy_cell_variables(const Ind fromIdx, const Ind toIdx) noexcept {
    mInt(toIdx) = mInt(fromIdx);
    mNum(toIdx) = mNum(fromIdx);
    for (SInd d = 0; d < nd; ++d) {
      mIntA(toIdx, d) = mIntA(fromIdx, d);
      mNumA(toIdx, d) = mNumA(fromIdx, d);
    }
  }
};

/// \brief Value type
template<SInd nd>
struct Value : container::sequential::ValueFacade<Container<nd>> {
  inline Int&      mInt()              noexcept { return mInt_;     }
  inline Num&      mNum()              noexcept { return mNum_;     }
  inline IntA<nd>& mIntA()             noexcept { return mIntA_;    }
  inline Int&      mIntA(const SInd d) noexcept { return mIntA_(d); }
  inline NumA<nd>& mNumA()             noexcept { return mNumA_;    }
  inline Num&      mNumA(const SInd d) noexcept { return mNumA_(d); }

  /// \brief Value Swap
  template<class Value1, class Value2>
  static inline void swap_values(Value1& lhs, Value2& rhs) noexcept {
    using std::swap;
    swap(lhs.mInt(), rhs.mInt());
    swap(lhs.mNum(), rhs.mNum());
    for (SInd d = 0; d < nd; ++d) {
      swap(lhs.mIntA(d), rhs.mIntA(d));
      swap(lhs.mNumA(d), rhs.mNumA(d));
    }
  }

  /// \brief Value Copy
  template<class Value1, class Value2>
  static inline void copy_values(Value1& from, Value2& to) noexcept {
    to.mInt() = from.mInt();
    to.mNum() = from.mNum();
    for (SInd d = 0; d < nd; ++d) {
      to.mIntA(d) = from.mIntA(d);
      to.mNumA(d) = from.mNumA(d);
    }
  }

 private:
  Int mInt_;
  Num mNum_;
  IntA<nd> mIntA_;
  NumA<nd> mNumA_;
};

/// \brief Reference type for the fixed-size sequential test container
template<SInd nd_>
struct Reference : container::sequential::ReferenceFacade<Container<nd_>> {
  static const SInd nd = nd_;
  using Base = container::sequential::ReferenceFacade<Container<nd>>;
  using Base::c;
  using Base::index;
  using Base::operator=;
  Reference<nd>& operator=(Reference rhs) { return Base::operator=(rhs); }

  inline Int& mInt()        noexcept { return mInt_; }
  inline Num& mNum()        noexcept { return mNum_; }
  inline Int& mIntA(SInd d) noexcept
  { return c() ? c()->mIntA(index(), d) : mIntA_(d); }
  inline Num& mNumA(SInd d) noexcept
  { return c() ? c()->mNumA(index(), d) : mNumA_(d); }

  Reference() = delete;
  Reference(Container<nd>* container, Ind idx)
      : Base(container, idx)
      , mInt_(c()->mInt(index()))
      , mNum_(c()->mNum(index()))
      , mIntA_(*static_cast<IntA<nd>*>(nullptr))
      , mNumA_(*static_cast<NumA<nd>*>(nullptr))
  {}
  Reference(Value<nd>& v)
    : Base()
    , mInt_(v.mInt())
    , mNum_(v.mNum())
    , mIntA_(v.mIntA())
    , mNumA_(v.mNumA())
  {}

 private:
  Int& mInt_;
  Num& mNum_;
  IntA<nd>& mIntA_;
  NumA<nd>& mNumA_;
};

}  // namespace fixed

////////////////////////////////////////////////////////////////////////////////
}  // namespace test
}  // namespace sequential
}  // namespace container
}  // namespace hom3
////////////////////////////////////////////////////////////////////////////////
#undef ENABLE_DBG_
////////////////////////////////////////////////////////////////////////////////
#endif
