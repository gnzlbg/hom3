#ifndef HOM3_CONTAINER_SEQUENTIAL_TEST_FIXEDCELL_HPP_
#define HOM3_CONTAINER_SEQUENTIAL_TEST_FIXEDCELL_HPP_
////////////////////////////////////////////////////////////////////////////////
// Includes:
#include "../../sequential.hpp"
/// Options:
#define ENABLE_DBG_ 0
#include "../../../misc/dbg.hpp"
////////////////////////////////////////////////////////////////////////////////

/// Forward declarations:
namespace test { namespace fixed {
template<SInd nd> struct Container;
template<SInd nd> struct Reference;
template<SInd nd> struct Value;
}} // test::fixed namespace

namespace container { namespace sequential {
/// Fixed Test Container trait:
template<SInd nd> struct traits<test::fixed::Container<nd>> {
  using container_type    = tag::fixed_nodes;
  using value_type = typename test::fixed::Value<nd>;
  using reference = typename test::fixed::Reference<nd>;
  using cell_index_type = Ind;
};
}} // container::sequential

namespace test { namespace fixed {

/// \brief Fixed-size sequential test container
template<SInd nd> struct Container : container::Sequential<Container<nd>> {

  /// Aliases:
  friend container::Sequential<Container<nd>>;
  template<template <SInd> class V_, SInd nd_ = 1>
  using M = container::Matrix<Container,container::matrix::tag::Cell,V_,Ind,SInd,nd_>;

  /// Data:
  M<IntM>    myInt;
  M<NumM>    myNum;
  M<IntM,nd> myIntArr;
  M<NumM,nd> myNumArr;

  /// Construction:
  Container(const Ind n)
      : container::Sequential< Container<nd> >(n,"fixed_container"),
      myInt(this,"myInt"), myNum(this,"myNum"), myIntArr(this,"myIntArr"),
      myNumArr(this,"myNumArr") {}

 private:

  inline void reset_cell(const Ind id) {
    myInt(id) = 0;
    myNum(id) = 0;
    for(SInd d = 0; d < nd; ++d) {
      myIntArr(id,d) = 0;
      myNumArr(id,d) = 0;
    }
  }

  inline void copy_cell_variables(const Ind fromId, const Ind toId) {
    myInt(toId) = myInt(fromId);
    myNum(toId) = myNum(fromId);
    for(SInd d = 0; d < nd; ++d) {
      myIntArr(toId,d) = myIntArr(fromId,d);
      myNumArr(toId,d) = myNumArr(fromId,d);
    }
  }

};


/// \brief Value type for the fixed-size sequential test container
template<SInd nd> struct Value : container::sequential::ValueFacade<Container<nd>> {
  inline Int& myInt() { return myInt_; }
  inline Num& myNum() { return myNum_; }
  inline IntA<nd>& myIntArr() { return myIntArr_; }
  inline Int& myIntArr(SInd d) { return myIntArr_(d); }
  inline NumA<nd>& myNumArr() { return myNumArr_; }
  inline Num& myNumArr(SInd d) { return myNumArr_(d); }

  template<class Value1, class Value2>
  static inline void swap_values(Value1& lhs, Value2& rhs) {
    std::swap(lhs.myInt(),rhs.myInt());
    std::swap(lhs.myNum(),rhs.myNum());
    for(SInd d = 0; d < nd; ++d) {
      std::swap(lhs.myIntArr(d),rhs.myIntArr(d));
      std::swap(lhs.myNumArr(d),rhs.myNumArr(d));
    }
  }

  template<class Value1, class Value2>
  static inline void copy_values(Value1& lhs, Value2& rhs) {
    lhs.myInt() = rhs.myInt();
    lhs.myNum() = rhs.myNum();
    for(SInd d = 0; d < nd; ++d) {
      lhs.myIntArr(d) = rhs.myIntArr(d);
      lhs.myNumArr(d) = rhs.myNumArr(d);
    }
  }
 private:
  Int myInt_;
  Num myNum_;
  IntA<nd> myIntArr_;
  NumA<nd> myNumArr_;
};


/// \brief Reference type for the fixed-size sequential test container
template<SInd nd> struct Reference : container::sequential::ReferenceFacade<Container<nd>> {
  using Base = container::sequential::ReferenceFacade<Container<nd>>;
  using Base::c;
  using Base::index;
  using Base::operator=;
  Reference<nd>& operator=(Reference rhs) { return Base::operator=(rhs); }

  inline Int& myInt() { return c()->myInt(index()); }
  inline Num& myNum() { return c()->myNum(index()); }
  inline Int& myIntArr(SInd d) { return c()->myIntArr(index(),d); }
  inline Num& myNumArr(SInd d) { return c()->myNumArr(index(),d); }

  Reference() : Base()  {}
  Reference(Container<nd>* c, Ind i) : Base(c,i) {}
};


}} // test::fixed namespace

#undef ENABLE_DBG_
////////////////////////////////////////////////////////////////////////////////
#endif
