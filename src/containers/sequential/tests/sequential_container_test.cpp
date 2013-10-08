////////////////////////////////////////////////////////////////////////////////
/// Library Includes:
#include <algorithm>
#include <vector>
#include <type_traits>
/// External Includes:
#include "gtest/gtest.h"
/// Includes:
#include "containers/sequential/tests/fixedcell.hpp"
#include "containers/sequential/tests/nodalcell.hpp"
/// Options:
#define ENABLE_DBG_ 0
#include "misc/dbg.hpp"
////////////////////////////////////////////////////////////////////////////////

using namespace hom3;

static const SInd nD = 2U;
using TestContainer2D = test::fixed::Container<2> ;
using TestNodalContainer2D = test::nodal::Container<2> ;

/// \test fixed_container constructor
///
/// Note: disabled until a decision is made about using/not
/// using exception handling
TEST(fixed_container_test, constructor) {
  // use error::terminate instead and EXPECT_DEATH
  // /// Constructor with zero size:
  // bool success = false;
  // try {
  //   TestContainer2D cells(0);
  // } catch (const std::exception& e) {
  //   success = true;
  // }
  //EXPECT_TRUE(success); ///< should throw a hom3 exception.

  /// TODO: fix this test
  // /// Constructor with zero size:
  // success = false;
  // try {
  //   TestContainer2D cells(-1); ///< dangerous: implicitly converts to unsigned!
  // } catch (const std::exception& e) {
  //   success = true;
  // }
  //EXPECT_TRUE(success); ///< should throw a bad_alloc.

  /// Empty constructor, copy constructor, and assignment
  /// should not compile.
  // TestContainer2D a;
  // TestContainer2D b(10);
  // TestContainer2D c(b);
  // TestContainer2D c = a;
}

template<class C> void initCellVariables2D(C& cells) {
  TRACE_IN_();
  for(Ind i = 0, e = cells.size(); i < e; ++i) {
    cells.myInt(i) = static_cast<Int>(i);
    cells.myNum(i) = static_cast<Num>(i) / 2;
    for(SInd d = 0; d < nD; ++d ) {
      cells.myIntArr(i,d) = static_cast<Int>(i+d);
      cells.myNumArr(i,d) = static_cast<Num>(i) / 2+d;
    }
  }
  TRACE_OUT();
}

template<class C> void checkCellVariables2D(C& cells, Ind fromId, Ind toId) {
  TRACE_IN_();
  for(Ind i = fromId, e = toId; i < e; ++i) {
    EXPECT_EQ(cells.myInt(i),static_cast<Int>(i));
    EXPECT_FLOAT_EQ(cells.myNum(i),static_cast<Num>(i) / 2);
    for(SInd d = 0; d < nD; ++d ) {
      EXPECT_EQ(cells.myIntArr(i,d),static_cast<Int>(i+d));
      EXPECT_FLOAT_EQ(cells.myNumArr(i,d),static_cast<Num>(i) / 2+d);
    }
  }
  TRACE_OUT();
}

void check_vars(TestContainer2D::reference cell) {
  EXPECT_EQ(cell.myInt(),static_cast<Int>(cell.index()));
  EXPECT_FLOAT_EQ(cell.myNum(),static_cast<Num>(cell.index())/2);
  for(SInd d = 0; d < nD; ++d ) {
    EXPECT_EQ(cell.myIntArr(d),static_cast<Int>(cell.index()+d));
    EXPECT_FLOAT_EQ(cell.myNumArr(d),static_cast<Num>(cell.index()) / 2+d);
  }
};


/// \test fixed_container iterators
TEST(fixed_container_test, iterators_DeathTest) {
  TestContainer2D cells(100);

  cells.push_cell(10);

  initCellVariables2D(cells);

  checkCellVariables2D(cells,0,cells.size());

  /// Read values (iterators,for_each):
  boost::for_each(cells, [](TestContainer2D::reference c){ check_vars(c); });

  for( auto i : cells ) {
    static_assert(std::is_same<decltype(i),TestContainer2D::reference>::value, "value type test"); 
  }

  /// Modify with random access:
  for(Ind i = 0, e = cells.size(); i < e; ++i) {
    cells[i].myInt() *= 2;
    cells[i].myNum() *= 3;
    for(SInd d = 0; d < nD; ++d ) {
      cells[i].myIntArr(d) *= 4;
      cells[i].myNumArr(d) *= 5;
    }
  }
  boost::for_each(cells,[](TestContainer2D::reference cell) {
      EXPECT_EQ(cell.myInt(),2*static_cast<Int>(cell.index()));
      EXPECT_FLOAT_EQ(cell.myNum(),3*static_cast<Num>(cell.index())/2);
      for(SInd d = 0; d < nD; ++d ) {
        EXPECT_EQ(cell.myIntArr(d),4*static_cast<Int>(cell.index()+d));
        EXPECT_FLOAT_EQ(cell.myNumArr(d),5*(static_cast<Num>(cell.index()) / 2+d));
      }
    });

  /// Access out of bounds elements: die in debug mode only
  EXPECT_DEBUG_DEATH(cells[-1].myNum() = -2.0, "[\\S\\s]+"); ///< implicit conversion to unsigned!
  EXPECT_DEBUG_DEATH(cells[cells.size()].myNum() = -2.0;, "[\\S\\s]+");
  EXPECT_DEBUG_DEATH(cells[cells.size()+1].myNum() = -2.0;, "[\\S\\s]+");
  EXPECT_DEBUG_DEATH(cells.myNum(-1) = -2.0, "[\\S\\s]+");
  EXPECT_DEBUG_DEATH(cells.myNum(cells.size()) = -2.0;, "[\\S\\s]+");
  EXPECT_DEBUG_DEATH(cells.myNum(cells.size()+1) = -2.0;, "[\\S\\s]+");
}

TEST(fixed_container_test, copy_cells) {
  TestContainer2D cells(100);

  cells.push_cell(10);

  initCellVariables2D(cells);

  /// Copy second half onto first half:
  for(Ind i = cells.size()/2, e = cells.size(), j = cells.size()/2; i < e; ++i) {
    cells.copy_cell(i,i-j);
  }

  // functions to check the solution:
  auto check_first_half = [&]() {
    for(Ind i = 0, e = cells.size()/2; i < e; ++i) {
      EXPECT_EQ(cells.myInt(i),cells.myInt(i+e));
      EXPECT_FLOAT_EQ(cells.myNum(i),cells.myNum(i+e));
      for(SInd d = 0; d < nD; ++d ) {
        EXPECT_EQ(cells.myIntArr(i,d),cells.myIntArr(i+e,d));
        EXPECT_FLOAT_EQ(cells.myNumArr(i,d),cells.myNumArr(i+e,d));
      }
    }
  };

  auto check_second_half = [&]() {
    checkCellVariables2D(cells,cells.size()/2,cells.size()); ///< Second half.
  };

  /// Test:
  check_first_half();
  check_second_half();

  /// Same test but with iterators:
  initCellVariables2D(cells);

  /// copy second half onto first half
  auto end = std::begin(cells) + cells.size()/2;
  std::copy(end,std::end(cells),std::begin(cells));

  /// Test:
  check_first_half();
  check_second_half();
}

TEST(fixed_container_test, erase_cells) {
  TestContainer2D cells(100);

  cells.push_cell(10);

  initCellVariables2D(cells);

  Ind oldSize = cells.size();

  // remove even elements
  auto next = boost::remove_if(cells,[](const TestContainer2D::reference c){
      return c.index() % 2 == 0;
    });

  EXPECT_EQ(next - std::begin(cells),5);

  for(Ind j = 0, ej = oldSize, i = 0; j < ej; ++j) {
    if(j % 2 == 0) { continue; }
    EXPECT_EQ(cells.myInt(i),static_cast<Int>(j));
    EXPECT_FLOAT_EQ(cells.myNum(i),static_cast<Num>(j) / 2);
    for(SInd d = 0; d < nD; ++d ) {
      EXPECT_EQ(cells.myIntArr(i,d),static_cast<Int>(j+d));
      EXPECT_FLOAT_EQ(cells.myNumArr(i,d),static_cast<Num>(j) / 2+d);
    }
    ++i;
  }
}

TEST(fixed_container_test, iterator_sort) {
  TestContainer2D cells(100);

  cells.push_cell(10);

  initCellVariables2D(cells);

  checkCellVariables2D(cells,0,cells.size());

  auto check_inversely_sorted = [&]() {
    for(Ind i = 0, e = 10; i < e; ++i) {
      EXPECT_EQ(cells.myInt(i),static_cast<Int>(9 - i));
      EXPECT_FLOAT_EQ(cells.myNum(i),static_cast<Num>(9 - i) / 2);
      for(SInd d = 0; d < nD; ++d ) {
        EXPECT_EQ(cells.myIntArr(i,d),static_cast<Int>(9 - i+d));
        EXPECT_FLOAT_EQ(cells.myNumArr(i,d),static_cast<Num>(9 - i) / 2+d);
      }
    }
  };

  boost::sort(cells,[](TestContainer2D::reference i, TestContainer2D::reference j){
      return i.myNum() > j.myNum();
    });

  check_inversely_sorted();
}


////////////////////////////////////////////////////////////////////////////////
/// Test Cells with variable number of nodes:
////////////////////////////////////////////////////////////////////////////////

template<class C> void plotCellNodes2D(C& cells) {
  TRACE_IN_();
  for(Ind cellId = 0, endCellId = cells.size(); cellId < endCellId; ++cellId) {
    DBG_ON("cId:",cellId,"[nB,nE]: [",cells.first_node(cellId),               \
        ",",cells.last_node(cellId),"]", "#nodes:",cells.node_size(cellId));
    plotCellNodeVariables2D(cells,cellId);
  }
  TRACE_OUT();
}

template<class C> void plotCellNodeVariables2D(C& cells,Ind cellId) {
  TRACE_IN((cellId));
  std::cerr << "# cId: " << cellId << " ";
  for(Ind nodeId = cells.first_node(cellId), endNodeId = cells.last_node(cellId);
      nodeId < endNodeId; ++nodeId)
    std::cerr << " | i:" << cells.myNodalInt(nodeId)
              <<  " n:" << cells.myNodalNum(nodeId);
  std::cerr << std::endl;
  TRACE_OUT();
}

template<class C> void initNodalVariables2D(C& cells) {
  TRACE_IN_();
  initCellVariables2D(cells);

  for(Ind cellId = 0, endCellId = cells.size(); cellId < endCellId; ++cellId) {
    for(Ind nodeId = cells.first_node(cellId),
         endNodeId = cells.last_node(cellId); nodeId < endNodeId; ++nodeId) {
      cells.myNodalInt(nodeId) = nodeId;
      cells.myNodalNum(nodeId) = 2*nodeId;
    }
  }
  TRACE_OUT();
}

template<class C> void checkNodalVariables2D(C& cells,C& cells_copy,Ind fromCellId,Ind toCellId) {
  TRACE_IN_();
  checkCellVariables2D(cells,fromCellId,toCellId);

  for(Ind cellId = fromCellId, endCellId = toCellId; cellId < endCellId; ++cellId) {
    for(Ind ni = cells.first_node(cellId), nei = cells.last_node(cellId),
            nc = cells_copy.first_node(cellId), nec = cells_copy.last_node(cellId);
            nc < nec; ++nc, ++ni) {
      EXPECT_EQ( (nei - ni), (nec - nc) );
      EXPECT_EQ(cells.myNodalInt(ni),cells_copy.myNodalInt(nc));
      EXPECT_FLOAT_EQ(cells.myNodalNum(ni),cells_copy.myNodalNum(nc));
    }
  }
  TRACE_OUT();
}

template<class C> void copy_last_half_to_beginning(C& cells, C& cells_copy) {
  TRACE_IN_();
  /// Copy cells:
  for(Ind i = cells.size()/2, e = cells.size(), j = cells.size()/2; i < e; ++i)
    cells.copy_cell(i,i-j);

  /// Test:
  for(Ind cellId = 0,
       endCellId = cells.size()/2; cellId < endCellId; ++cellId) { ///< First half.
    EXPECT_EQ(cells.myInt(cellId),cells_copy.myInt(cellId+endCellId));
    EXPECT_FLOAT_EQ(cells.myNum(cellId),cells_copy.myNum(cellId+endCellId));
    for(SInd d = 0; d < nD; ++d ) {
      EXPECT_EQ(cells.myIntArr(cellId,d),cells_copy.myIntArr(cellId+endCellId,d));
      EXPECT_FLOAT_EQ(cells.myNumArr(cellId,d),cells_copy.myNumArr(cellId+endCellId,d));
    }
    for(Ind ni = cells.first_node(cellId), nei = cells.last_node(cellId),
            nc = cells_copy.first_node(cellId+endCellId), nec = cells_copy.last_node(cellId+endCellId);
        nc < nec; ++nc, ++ni) {
      EXPECT_TRUE(ni < nei);
      EXPECT_EQ(cells.myNodalInt(ni),cells_copy.myNodalInt(nc));
      EXPECT_FLOAT_EQ(cells.myNodalNum(ni),cells_copy.myNodalNum(nc));
    }
  }

  checkNodalVariables2D(cells,cells_copy,cells.size()/2,cells.size()); ///< Second half.
  TRACE_OUT();
}

TestNodalContainer2D create_linear_increasing_nodal_2d_container(Ind noCells,Ind maxNoCells) {

  Ind linearMaxNoNodes = noCells*(noCells+1)/2;

  TestNodalContainer2D cells(maxNoCells,2*linearMaxNoNodes);
  for(Ind cellId = 0; cellId < noCells; ++cellId)
    cells.push_cell(1,cellId+1);
  initNodalVariables2D(cells);

  EXPECT_EQ(cells.size(),noCells);
  EXPECT_EQ(cells.node_size(),linearMaxNoNodes);
  return cells;
}

TestNodalContainer2D create_linear_decreasing_nodal_2d_container(Ind noCells,Ind maxNoCells) {

  Ind linearMaxNoNodes = noCells*(noCells+1)/2;

  TestNodalContainer2D cells(maxNoCells,2*linearMaxNoNodes);
  for(Ind cellId = noCells; cellId > 0; --cellId)
    cells.push_cell(1,cellId);
  initNodalVariables2D(cells);

  EXPECT_EQ(cells.size(),noCells);
  EXPECT_EQ(cells.node_size(),linearMaxNoNodes);
  return cells;
}

template<class C> Ind erase_even_cells_and_verify(C& cells, C& cells_copy) {
  TRACE_IN_();
  /// Erase even cells:
  auto oldSize = cells.size();
  auto pred = [](TestNodalContainer2D::reference c){ return c.index() % 2 == 0; };
  auto c = container::sequential::algorithm::erase_remove_if(cells,pred);
  // cells.pop_cell(std::end(cells) - next);
  //auto next = cells.erase_remove_if(pred);
  /// Check:
  for(TestNodalContainer2D::cell_size_type j = 0, ej = oldSize, i = 0; j < ej; ++j) {
    if(j % 2 == 0) continue;
    /// cell values:
    EXPECT_EQ(cells.myInt(i),cells_copy.myInt(j));
    EXPECT_FLOAT_EQ(cells.myNum(i),cells_copy.myNum(j));
    for(SInd d = 0; d < nD; ++d ) {
      EXPECT_EQ(cells.myIntArr(i,d),cells_copy.myIntArr(j,d));
      EXPECT_FLOAT_EQ(cells.myNumArr(i,d),cells_copy.myNumArr(j,d));
    }
    /// nodal values:
    for(Ind ni = cells.first_node(i), nei = cells.last_node(i),
            nc = cells_copy.first_node(j), nec = cells_copy.last_node(j);
                 nc < nec; ++nc, ++ni) {
      EXPECT_TRUE(ni < nei);
      EXPECT_EQ(cells.myNodalInt(ni),cells_copy.myNodalInt(nc));
      EXPECT_FLOAT_EQ(cells.myNodalNum(ni),cells_copy.myNodalNum(nc));
    }
    ++i;
  }
  TRACE_OUT();
  return c.size();
}

/// Nodal container where all cells have the same #of nodes:
TEST(variable_container_test, constant_copy) {

  Ind noCells = 10;
  Ind maxNoCells = 10, maxNoNodes = 100;

  /// Constant container:
  Ind nodesPerCell = 10;

  TestNodalContainer2D cells(maxNoCells,maxNoNodes);
  cells.push_cell(noCells,nodesPerCell);
  initNodalVariables2D(cells);

  TestNodalContainer2D cells_copy(cells);

  checkNodalVariables2D(cells,cells_copy,0,cells.size());

  copy_last_half_to_beginning(cells,cells_copy);
}

TEST(variable_container_test, constant_erase) {
  Ind noCells = 10;
  Ind maxNoCells = 10, maxNoNodes = 100;

  /// Constant container:
  Ind nodesPerCell = 10;

  /// Create container and copy:
  TestNodalContainer2D cells(maxNoCells,maxNoNodes);
  cells.push_cell(noCells,nodesPerCell);
  initNodalVariables2D(cells);
  TestNodalContainer2D cells_copy(cells);
  checkNodalVariables2D(cells,cells_copy,0,cells.size());

  /// Erase even cells:
  auto size = erase_even_cells_and_verify(cells,cells_copy);

  /// Check size:
  EXPECT_EQ(size,maxNoCells/2);

}

/// Nodal container where the #of nodes in the cells is incremented linearly
TEST(variable_container_test, linear_incr_copy) {
  Ind noCells = 10;
  Ind maxNoCells = 10; // maxNoNodes = 100;

  TestNodalContainer2D cells(create_linear_increasing_nodal_2d_container(noCells,maxNoCells));
  TestNodalContainer2D cells_copy(cells);
  checkNodalVariables2D(cells,cells_copy,0,cells.size());

  copy_last_half_to_beginning(cells,cells_copy);
}

TEST(variable_container_test, linear_incr_erase) {

  Ind noCells = 10;
  Ind maxNoCells = 10; // maxNoNodes = 100;

  TestNodalContainer2D cells(create_linear_increasing_nodal_2d_container(noCells,maxNoCells));
  TestNodalContainer2D cells_copy(cells);
  checkNodalVariables2D(cells,cells_copy,0,cells.size());

  auto size = erase_even_cells_and_verify(cells,cells_copy);

  /// Check size:
  EXPECT_EQ(size,maxNoCells/2);
  EXPECT_EQ(cells.size(),static_cast<Ind>(maxNoCells/2));
}

/// Nodal container where the #of nodes in the cells is decremented linearly
TEST(variable_container_test, linear_decr_copy) {

  Ind noCells = 10;
  Ind maxNoCells = 10; // maxNoNodes = 100;

  TestNodalContainer2D cells(create_linear_decreasing_nodal_2d_container(noCells,maxNoCells));
  TestNodalContainer2D cells_copy(cells);
  checkNodalVariables2D(cells,cells_copy,0,cells.size());

  copy_last_half_to_beginning(cells,cells_copy);
}

TEST(variable_container_test, linear_decr_erase) {

  Ind noCells = 10;
  Ind maxNoCells = 10; // maxNoNodes = 100;

  TestNodalContainer2D cells(create_linear_decreasing_nodal_2d_container(noCells,maxNoCells));
  TestNodalContainer2D cells_copy(cells);
  checkNodalVariables2D(cells,cells_copy,0,cells.size());

  auto size = erase_even_cells_and_verify(cells,cells_copy);

  /// Check size:
  EXPECT_EQ(size,(maxNoCells/2));
  EXPECT_EQ(cells.size(),static_cast<Ind>(maxNoCells/2));
}


TEST(variable_container_test, single_index_storage) {

  Ind maxNoCells = 5, maxNoNodes = 20;

  TestNodalContainer2D cells(maxNoCells,maxNoNodes);
  using nst = TestNodalContainer2D::node_size_type;
  cells.push_cell(1,1);  cells.push_cell(1,2);  cells.push_cell(1,3);

  //plotCellNodes2D(cells);
  DBG("");
  EXPECT_EQ(cells.node_size(0), nst{1});
  EXPECT_EQ(cells.node_size(1), nst{2});
  EXPECT_EQ(cells.node_size(2), nst{3});

  DBG("cells.erase_remove_if([&cells](Ind i)->bool{ return i == 2; });");
  cells.erase_remove_if([&cells](Ind i)->bool{ return i == 2; });
  //plotCellNodes2D(cells);
  EXPECT_EQ(cells.node_size(0), nst{1});
  EXPECT_EQ(cells.node_size(1), nst{2});
  DBG("");

  DBG("cells.push_back(1,5);");
  cells.push_cell(1,5);
  //plotCellNodes2D(cells);
  EXPECT_EQ(cells.node_size(0), nst{1});
  EXPECT_EQ(cells.node_size(1), nst{2});
  EXPECT_EQ(cells.node_size(2), nst{5});
  DBG("");

  DBG("cells.erase_remove_if([&cells](Ind i)->bool{ return i == 0; });");
  cells.erase_remove_if([&cells](Ind i)->bool{ return i == 0; });
  //plotCellNodes2D(cells);
  DBG("");

  DBG("cells.push_back(1,3);");
  cells.push_cell(1,3);
  //plotCellNodes2D(cells);
  DBG("");

  DBG("cells.erase_remove_if([&cells](Ind i)->bool{ return i == 2; });");
  cells.erase_remove_if([&cells](Ind i)->bool{ return i == 2; });
  //plotCellNodes2D(cells);
  DBG("");

  DBG("cells.push_back(1,4);");
  cells.push_cell(1,4);
  //plotCellNodes2D(cells);
  DBG("");

}

////////////////////////////////////////////////////////////////////////////////
#undef ENABLE_DBG_
////////////////////////////////////////////////////////////////////////////////
