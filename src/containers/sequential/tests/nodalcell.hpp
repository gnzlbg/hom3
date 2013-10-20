#ifndef HOM3_CONTAINER_SEQUENTIAL_TEST_NODALCELL_HPP_
#define HOM3_CONTAINER_SEQUENTIAL_TEST_NODALCELL_HPP_
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

/// Nodal-container Forward declarations:
namespace nodal {
template<SInd nd> struct Container;
template<SInd nd> struct Reference;
template<SInd nd> struct Value;
}  // namespace nodal

}  // namespace test

/// \brief Nodal-container trait
template<SInd nd> struct traits<test::nodal::Container<nd>> {
  using container_type  = tag::variable_nodes;
  using value_type      = typename test::nodal::Value<nd>;
  using reference       = typename test::nodal::Reference<nd>;
  using cell_index_type = Ind;
  using node_index_type = Ind;
};

namespace test { namespace nodal {

/// \brief Test sequential container with variable #of nodes per element
template<SInd nd_> struct Container : container::Sequential<Container<nd_>> {
  static const SInd nd = nd_;
  friend container::Sequential<Container<nd>>;

  /// Aliases:
  /// Cell Matrix
  template<template <SInd> class V_, SInd nd__ = 1> using CM = container::Matrix
  <Container, container::matrix::tag::Cell, V_, Ind, SInd, nd__>;
  /// Nodal Matrix
  template<template <SInd> class V_, SInd nd__ = 1> using NM = container::Matrix
  <Container, container::matrix::tag::Node, V_, Ind, SInd, nd__>;

  /// Construction:
  Container(const Ind ne, const Ind nn) noexcept
    : container::Sequential<Container<nd>>(ne, nn, "nodal_container")
    ,  mInt(this, "mInt")
    ,  mNum(this, "mNum")
    , mIntA(this, "mIntA")
    , mNumA(this, "mNumA")
    , mIntN(this, "mIntN")
    , mNumN(this, "mNumN")
  {}

  /// Data:
  CM<IntM>     mInt;
  CM<NumM>     mNum;
  CM<IntM, nd> mIntA;
  CM<NumM, nd> mNumA;
  NM<IntM>     mIntN;
  NM<NumM>     mNumN;

  /// \name Extra functionality, not required as long as you don't use it!
  ///@{

  /// \brief Copy/Assing/Move containers:
  Container(const Container<nd>& other)
    : container::Sequential<Container<nd>>(other)
    , mInt(other.mInt)
    , mNum(other.mNum)
    , mIntA(other.mIntA)
    , mNumA(other.mNumA)
    , mIntN(other.mIntN)
    , mNumN(other.mNumN) {
    TRACE_IN_();

    this->initialize(mInt, mNum, mIntA, mNumA, mIntN, mNumN);

    ASSERT(this->size()      == other.size()        , "Wrong size!");
    ASSERT(this->node_size() == other.node_size()   , "Wrong node size!");
    ASSERT(mInt().size()     == other.mInt().size() , "Wrong size!");
    ASSERT(mNum().size()     == other.mNum().size() , "Wrong size!");
    ASSERT(mIntA().size()    == other.mIntA().size(), "Wrong size!");
    ASSERT(mNumA().size()    == other.mNumA().size(), "Wrong size!");
    ASSERT(mIntN().size()    == other.mIntN().size(), "Wrong size!");
    ASSERT(mNumN().size()    == other.mNumN().size(), "Wrong size!");

    TRACE_OUT();
  }
  Container& operator=(Container<nd> other) noexcept {
    TRACE_IN_();
    swap(*this, other);
    this->initialize(mInt, mNum, mIntA, mNumA, mIntN, mNumN);
    TRACE_OUT();
    return *this;
  }
  Container(Container<nd>&& other) noexcept : Container(0, 0) {
    TRACE_IN_();
    swap(*this, other);
    this->initialize(mInt, mNum, mIntA, mNumA, mIntN, mNumN);
    TRACE_OUT();
  }
  ///@}

 private:
  inline void reset_cell(const Ind cIdx) noexcept {
    TRACE_IN((cIdx));
    mInt(cIdx) = 0;
    mNum(cIdx) = 0;
    for (SInd d = 0; d < nd; ++d) {
      mIntA(cIdx, d) = 0;
      mNumA(cIdx, d) = 0;
    }
    TRACE_OUT();
  }

  inline
  void copy_cell_variables(const Ind fromCIdx, const Ind toCIdx) noexcept {
    TRACE_IN((fromCIdx)(toCIdx));
    mInt(toCIdx) = mInt(fromCIdx);
    mNum(toCIdx) = mNum(fromCIdx);
    for (SInd d = 0; d < nd; ++d) {
      mIntA(toCIdx, d) = mIntA(fromCIdx, d);
      mNumA(toCIdx, d) = mNumA(fromCIdx, d);
    }
    TRACE_OUT();
  }

  inline
  void swap_cell_variables(const Ind fromCIdx, const Ind toCIdx) noexcept {
    TRACE_IN((fromCIdx)(toCIdx));
    using std::swap;
    swap(mInt(toCIdx), mInt(fromCIdx));
    swap(mNum(toCIdx), mNum(fromCIdx));
    for (SInd d = 0; d < nd; ++d) {
      swap(mIntA(toCIdx, d), mIntA(fromCIdx, d));
      swap(mNumA(toCIdx, d), mNumA(fromCIdx, d));
    }
    TRACE_OUT();
  }

  inline
  void copy_node_variables(const Ind fromNIdx, const Ind toNIdx) noexcept {
    mIntN(toNIdx) = mIntN(fromNIdx);
    mNumN(toNIdx) = mNumN(fromNIdx);
  }

  inline
  void swap_node_variables(const Ind fromNIdx, const Ind toNIdx) noexcept {
    using std::swap;
    swap(mIntN(toNIdx), mIntN(fromNIdx));
    swap(mNumN(toNIdx), mNumN(fromNIdx));
  }

  /// \name Extra functionality, not required as long as you don't use it!
  ///@{

  /// Required for Copy/Assing/Move containers only:
  friend
  void swap_containers(Container<nd>& first, Container<nd>& second) noexcept {
    TRACE_IN_();
    using std::swap;
    swap(first.mInt, second.mInt);
    swap(first.mNum, second.mNum);
    swap(first.mIntA, second.mIntA);
    swap(first.mNumA, second.mNumA);
    swap(first.mIntN, second.mIntN);
    swap(first.mNumN, second.mNumN);
    TRACE_OUT();
  };
  ///@}
};

/// \brief Value-type for variable node container
template<SInd nd>
struct Value : container::sequential::ValueFacade<Container<nd>> {};

/// Reference Type for variable node container
template<SInd nd>
struct Reference : container::sequential::ReferenceFacade<Container<nd>> {
  using Base = container::sequential::ReferenceFacade<Container<nd>>;
  using Base::c;
  using Base::index;
  using Base::operator=;
  Reference<nd>& operator=(Reference rhs) { return Base::operator=(rhs); }

  Reference() : Base()  {}
  Reference(Container<nd>* c, Ind i) : Base(c, i) {}
};

}  // namespace nodal

////////////////////////////////////////////////////////////////////////////////
}  // namespace test
}  // namespace sequential
}  // namespace container
}  // namespace hom3
////////////////////////////////////////////////////////////////////////////////
#undef ENABLE_DBG_
////////////////////////////////////////////////////////////////////////////////
#endif
