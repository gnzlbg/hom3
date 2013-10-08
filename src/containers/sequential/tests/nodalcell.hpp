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
namespace hom3 {
////////////////////////////////////////////////////////////////////////////////

/// Forward declarations:
namespace test { namespace nodal {
template<SInd nd> struct Container;
template<SInd nd> struct Reference;
template<SInd nd> struct Value;
}} // test::nodal namespace

namespace container { namespace sequential {
/// \brief Specialize trait for nodal test container
template<SInd nd> struct traits<test::nodal::Container<nd>> {
  using container_type    = tag::variable_nodes;
  using value_type        = typename test::nodal::Value<nd>;
  using reference         = typename test::nodal::Reference<nd>;
  using cell_index_type   = Ind;
  using node_index_type   = Ind;
};
}} // container::sequential namespace

namespace test { namespace nodal {

/// \brief Test sequential container with variable #of nodes per element
template<SInd nd> struct Container : container::Sequential<Container<nd>> {

  friend container::Sequential<Container<nd>>;

  /// Aliases:

  /// Cell Matrix
  template<template <SInd> class V_, SInd nd_ = 1>
      using CM = container::Matrix<Container,container::matrix::tag::Cell,V_,Ind,SInd,nd_>;

  /// Nodal Matrix
  template<template <SInd> class V_, SInd nd_ = 1>
      using NM = container::Matrix<Container,container::matrix::tag::Node,V_,Ind,SInd,nd_>;

  /// Data:
  CM<IntM>    myInt;
  CM<NumM>    myNum;
  CM<IntM,nd> myIntArr;
  CM<NumM,nd> myNumArr;
  NM<IntM>    myNodalInt;
  NM<NumM>    myNodalNum;

  /// Construction:
  Container(const Ind ne, const Ind nn)
      : container::Sequential<Container<nd>>(ne,nn,"nodal_container"),
      myInt(this,"myInt"), myNum(this,"myNum"), myIntArr(this,"myIntArr"),
      myNumArr(this,"myNumArr"), myNodalInt(this,"myNodalInt"),
      myNodalNum(this,"myNodalNum") {}

  /// \name Extra functionality, not required as long as you don't use it!
  ///@{

  /// \brief Copy/Assing/Move containers:
  Container(const Container<nd>& other)
      : container::Sequential<Container<nd>>(other),
        myInt(other.myInt), myNum(other.myNum),
        myIntArr(other.myIntArr), myNumArr(other.myNumArr),
        myNodalInt(other.myNodalInt), myNodalNum(other.myNodalNum) {
    TRACE_IN_();

    this->initialize(myInt,myNum,myIntArr,myNumArr,myNodalInt,myNodalNum);

    ASSERT(this->size()        == other.size()              ,"Wrong size!");
    ASSERT(this->node_size()   == other.node_size()         ,"Wrong node size!");
    ASSERT(myInt().size()      == other.myInt().size()      ,"Wrong size!");
    ASSERT(myNum().size()      == other.myNum().size()      ,"Wrong size!");
    ASSERT(myIntArr().size()   == other.myIntArr().size()   ,"Wrong size!");
    ASSERT(myNumArr().size()   == other.myNumArr().size()   ,"Wrong size!");
    ASSERT(myNodalInt().size() == other.myNodalInt().size() ,"Wrong size!");
    ASSERT(myNodalNum().size() == other.myNodalNum().size() ,"Wrong size!");

    TRACE_OUT();
  }
  Container& operator=(Container<nd> other) {
    TRACE_IN_();
    swap(*this, other);
    this->initialize(myInt,myNum,myIntArr,myNumArr,myNodalInt,myNodalNum);
    TRACE_OUT();
    return *this;
  }
  Container(Container<nd>&& other) : Container(0,0) {
    TRACE_IN_();
    swap(*this, other);
    this->initialize(myInt,myNum,myIntArr,myNumArr,myNodalInt,myNodalNum);
    TRACE_OUT();
  }
  ///@}

 private:

  inline void reset_cell(const Ind id) {
    TRACE_IN((id));
    myInt(id) = 0;
    myNum(id) = 0;
    for(SInd d = 0; d < nd; ++d) {
      myIntArr(id,d) = 0;
      myNumArr(id,d) = 0;
    }
    TRACE_OUT();
  }

  inline void copy_cell_variables(const Ind fromCellId, const Ind toCellId) {
    TRACE_IN((fromCellId)(toCellId));
    myInt(toCellId) = myInt(fromCellId);
    myNum(toCellId) = myNum(fromCellId);
    for(SInd d = 0; d < nd; ++d) {
      myIntArr(toCellId,d) = myIntArr(fromCellId,d);
      myNumArr(toCellId,d) = myNumArr(fromCellId,d);
    }
    TRACE_OUT();
  }

  inline void swap_cell_variables(const Ind fromCellId, const Ind toCellId) {
    TRACE_IN((fromCellId)(toCellId));
    using std::swap;
    swap(myInt(toCellId),myInt(fromCellId));
    swap(myNum(toCellId),myNum(fromCellId));
    for(SInd d = 0; d < nd; ++d) {
      swap(myIntArr(toCellId,d),myIntArr(fromCellId,d));
      swap(myNumArr(toCellId,d),myNumArr(fromCellId,d));
    }
    TRACE_OUT();
  }

  inline void copy_node_variables(const Ind fromNodeId, const Ind toNodeId) {
    myNodalInt(toNodeId) = myNodalInt(fromNodeId);
    myNodalNum(toNodeId) = myNodalNum(fromNodeId);
  }

  inline void swap_node_variables(const Ind fromNodeId, const Ind toNodeId) {
    using std::swap;
    swap(myNodalInt(toNodeId),myNodalInt(fromNodeId));
    swap(myNodalNum(toNodeId),myNodalNum(fromNodeId));
  }


  /// \name Extra functionality, not required as long as you don't use it!
  ///@{

  /// Required for Copy/Assing/Move containers only:
  friend void swap_containers (Container<nd>& first, Container<nd>& second) {
    TRACE_IN_();
    using std::swap;
    swap(first.myInt      , second.myInt      );
    swap(first.myNum      , second.myNum      );
    swap(first.myIntArr   , second.myIntArr   );
    swap(first.myNumArr   , second.myNumArr   );
    swap(first.myNodalInt , second.myNodalInt );
    swap(first.myNodalNum , second.myNodalNum );
    TRACE_OUT();
  };
  ///@}

};


/// \brief Value-type for variable node container
template<SInd nd> struct Value : container::sequential::ValueFacade<Container<nd>> {};

/// Reference Type for variable node container
template<SInd nd> struct Reference : container::sequential::ReferenceFacade<Container<nd>> {
  using Base = container::sequential::ReferenceFacade<Container<nd>>;
  using Base::c;
  using Base::index;
  using Base::operator=;
  Reference<nd>& operator=(Reference rhs) { return Base::operator=(rhs); }

  Reference() : Base()  {}
  Reference(Container<nd>* c, Ind i) : Base(c,i) {}
};

}}} // hom3::test::nodal namespace

#undef ENABLE_DBG_
////////////////////////////////////////////////////////////////////////////////
#endif
