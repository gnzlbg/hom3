/// \file \brief Tests for container::sequential
/// Library Includes:
#include <algorithm>
#include <vector>
#include <type_traits>
/// External Includes:
#include "misc/test.hpp"
/// Includes:
#include "containers/sequential/tests/helpers.hpp"
/// Options:
#define ENABLE_DBG_ 0
#include "misc/dbg.hpp"
////////////////////////////////////////////////////////////////////////////////
using namespace hom3;
using namespace container::sequential::test;

/// \test Test constructors
TEST(container_tests, constructors) {
  auto a = FC2D{10};
  auto b = NC2D{10, 50};

  /// Can't construct contaienrs of empty or negative sizes
  EXPECT_DEATH_(FC2D{0});
  EXPECT_DEATH_(FC2D{-1});

  EXPECT_DEATH_((NC2D{0, 0}));
  EXPECT_DEATH_((NC2D{-1, -1}));

  /// Nodal containers must have more nodes than elements
  EXPECT_DEATH_((NC2D{5, 2}));
}

/// \test Iterators, values, and references
TEST(container_tests, values_and_references) {
  /// Create a container with just one element
  auto cells = FC2D{1};
  cells.push_cell();
  init_cell<2>(cells[0], 1);

  /// Get a reference to that element by dereferencing an iterator
  auto ref = *std::begin(cells);
  static_assert(std::is_same<decltype(ref), FC2DRef>::value,
                "Dereferencing an iterator must return a reference!");
  check_cell_equality<2>(ref, *std::begin(cells));
  EXPECT_EQ(ref.mInt(), cells.mInt(0));
  EXPECT_NUM_EQ(ref.mNum(), cells.mNum(0));

  /// Modify the container element through the reference
  ref.mInt() += 1;
  ref.mNum() += 1.;
      EXPECT_EQ(ref.mInt(), 2);
  EXPECT_NUM_EQ(ref.mNum(), 3.);
      EXPECT_EQ(ref.mInt(), cells.mInt(0));
  EXPECT_NUM_EQ(ref.mNum(), cells.mNum(0));
  for (SInd d = 0; d < 2; ++d) {
    ref.mIntA(d) += 1;
    ref.mNumA(d) += 1.;
        EXPECT_EQ(ref.mIntA(d), cells.mIntA(0, d));
    EXPECT_NUM_EQ(ref.mNumA(d), cells.mNumA(0, d));
  }

  /// Make a copy of the container element reference into a value
  FC2DVal val = ref;
      EXPECT_EQ(val.mInt(), cells.mInt(0));
  EXPECT_NUM_EQ(val.mNum(), cells.mNum(0));
  for (SInd d = 0; d < 2; ++d) {
        EXPECT_EQ(val.mIntA(d), cells.mIntA(0, d));
    EXPECT_NUM_EQ(val.mNumA(d), cells.mNumA(0, d));
  }
  /// Modify the value does not modify the container element
  val.mInt() += 1;
  val.mNum() += 1.;
      EXPECT_EQ(val.mInt(), 3);
  EXPECT_NUM_EQ(val.mNum(), 4.);
      EXPECT_EQ(val.mInt(), cells.mInt(0) + 1);
  EXPECT_NUM_EQ(val.mNum(), cells.mNum(0) + 1.);
  for (SInd d = 0; d < 2; ++d) {
    val.mIntA(d) += 1;
    val.mNumA(d) += 1.;
        EXPECT_EQ(val.mIntA(d), cells.mIntA(0, d) + 1);
    EXPECT_NUM_EQ(val.mNumA(d), cells.mNumA(0, d) + 1.);
  }

  /// Asign modified value to reference modifies the container
  ref = val;
      EXPECT_EQ(val.mInt(), cells.mInt(0));
  EXPECT_NUM_EQ(val.mNum(), cells.mNum(0));
  for (SInd d = 0; d < 2; ++d) {
        EXPECT_EQ(val.mIntA(d), cells.mIntA(0, d));
    EXPECT_NUM_EQ(val.mNumA(d), cells.mNumA(0, d));
  }

  /// Modify the value via a reference to the value
  FC2DRef refToVal = val;
  refToVal.mInt() += 1;
  refToVal.mNum() += 1.;
  for (SInd d = 0; d < 2; ++d) {
    val.mIntA(d) += 1;
    val.mNumA(d) += 1.;
  }
      EXPECT_EQ(refToVal.mInt(), 4);
  EXPECT_NUM_EQ(refToVal.mNum(), 5.);
      EXPECT_EQ(val.mInt(), 4);
  EXPECT_NUM_EQ(val.mNum(), 5.);
      EXPECT_EQ(val.mInt(), cells.mInt(0) + 1);
  EXPECT_NUM_EQ(val.mNum(), cells.mNum(0) + 1.);
  for (SInd d = 0; d < 2; ++d) {
        EXPECT_EQ(refToVal.mIntA(d), cells.mIntA(0, d) + 1);
    EXPECT_NUM_EQ(refToVal.mNumA(d), cells.mNumA(0, d) + 1.);
        EXPECT_EQ(val.mIntA(d), cells.mIntA(0, d) + 1);
    EXPECT_NUM_EQ(val.mNumA(d), cells.mNumA(0, d) + 1.);
  }

  /// Can assign a reference to value to a new value
  FC2DVal newVal = refToVal;
      EXPECT_EQ(val.mInt(), newVal.mInt());
  EXPECT_NUM_EQ(val.mNum(), newVal.mNum());

  /// Can assign a reference to value to a reference into a container element,
  /// which modifies the container element
  ref = refToVal;
      EXPECT_EQ(val.mInt(), cells.mInt(0));
  EXPECT_NUM_EQ(val.mNum(), cells.mNum(0));
  for (SInd d = 0; d < 2; ++d) {
        EXPECT_EQ(val.mIntA(d), cells.mIntA(0, d));
    EXPECT_NUM_EQ(val.mNumA(d), cells.mNumA(0, d));
  }

  /// Assignment between different containers
  auto cells2 = FC2D{1};
  cells2.push_cell();
  cells2.mInt(0) = 1;
  cells2.mNum(0) = 2;
  *std::begin(cells) = *std::begin(cells2);
      EXPECT_EQ(cells.mInt(0), cells2.mInt(0));
  EXPECT_NUM_EQ(cells.mNum(0), cells2.mNum(0));

  *std::begin(cells) = val;
      EXPECT_EQ(cells.mInt(0), 4);
  EXPECT_NUM_EQ(cells.mNum(0), 5.);

  /// References from two values:
  FC2DVal val1, val2;
  FC2DRef val1r = val1, val2r = val2;
  val1.mInt() = 1;
  val1.mNum() = 2.;
  val2.mInt() = 3;
  val2.mNum() = 4.;
  val2r = val1r;
      EXPECT_EQ(val2.mInt(), val1.mInt());
  EXPECT_NUM_EQ(val2.mNum(), val2.mNum());
}

/// \test fixed_container iterators
TEST(fixed_container_test, iterators_DeathTest) {
  auto cells = FC2D{100};

  /// Initialize a container with 10 cells:
  cells.push_cell(10);
  init_variables(cells);
  check_variables(cells);

  /// Read values (iterators,for_each):
  algorithm::for_each(cells, [](FC2DRef c) { check_vars(c); });

  /// Modify with random access:
  for (Ind i = 0, e = cells.size(); i < e; ++i) {
    cells[i].mInt() *= 2;
    cells[i].mNum() *= 3;
    for (SInd d = 0; d < FC2D::nd; ++d) {
      cells[i].mIntA(d) *= 4;
      cells[i].mNumA(d) *= 5;
    }
  }

  /// Recheck using references
  algorithm::for_each(cells, [](FC2DRef cell) {
    const auto i = cell.index();
        EXPECT_EQ(cell.mInt(), 2 * static_cast<Int>(i));
    EXPECT_NUM_EQ(cell.mNum(), 3 * static_cast<Num>(i) / 2);
    for (SInd d = 0; d < FC2D::nd; ++d) {
          EXPECT_EQ(cell.mIntA(d), 4 * static_cast<Int>(i + d));
      EXPECT_NUM_EQ(cell.mNumA(d), 5 * (static_cast<Num>(i) / 2 + d));
    }
  });

  /// Access out of bounds elements: die in debug mode only
  EXPECT_DEBUG_DEATH_(cells[-1].mNum()               = -2.0);
  EXPECT_DEBUG_DEATH_(cells[cells.size()].mNum()     = -2.0);
  EXPECT_DEBUG_DEATH_(cells[cells.size() + 1].mNum() = -2.0);
  EXPECT_DEBUG_DEATH_(cells.mNum(-1)                 = -2.0);
  EXPECT_DEBUG_DEATH_(cells.mNum(cells.size())       = -2.0);
  EXPECT_DEBUG_DEATH_(cells.mNum(cells.size() + 1)   = -2.0);
}

TEST(fixed_container_test, copy_cells) {
  FC2D cells(100);
  cells.push_cell(10);
  init_variables(cells);

  /// Copy second half onto first half:
  for (Ind i = cells.size() / 2, e = cells.size(),
           j = cells.size() / 2; i < e; ++i) {
    cells.copy_cell(i, i - j);
  }

  // functions to check the solution:
  auto check_first_half = [&]() {
    for (Ind i = 0, e = cells.size() / 2; i < e; ++i) {
      EXPECT_EQ(cells.mInt(i), cells.mInt(i + e));
      EXPECT_NUM_EQ(cells.mNum(i), cells.mNum(i + e));
      for (SInd d = 0; d < FC2D::nd; ++d) {
        EXPECT_EQ(cells.mIntA(i, d), cells.mIntA(i + e, d));
        EXPECT_NUM_EQ(cells.mNumA(i, d), cells.mNumA(i + e, d));
      }
    }
  };

  auto check_second_half = [&]() {  ///< Second half.
    check_variables(cells, cells.size() / 2, cells.size());
  };

  /// Test:
  check_first_half();
  check_second_half();

  /// Same test but with iterators:
  init_variables(cells);

  /// copy second half onto first half
  auto end = std::begin(cells) + cells.size() / 2;
  std::copy(end, std::end(cells), std::begin(cells));

  /// Test:
  check_first_half();
  check_second_half();
}

TEST(fixed_container_test, erase_cells) {
  FC2D cells(100);
  cells.push_cell(10);
  init_variables(cells);

  Ind oldSize = cells.size();

  // remove even elements
  auto next = algorithm::remove_if(cells, [](const FC2DRef c) {
    return c.index() % 2 == 0;
  });

  EXPECT_EQ(next - std::begin(cells), 5);

  for (Ind j = 0, ej = oldSize, i = 0; j < ej; ++j) {
    if (j % 2 == 0) { continue; }
        EXPECT_EQ(cells.mInt(i), static_cast<Int>(j));
    EXPECT_NUM_EQ(cells.mNum(i), static_cast<Num>(j) / 2);
    for (SInd d = 0; d < FC2D::nd; ++d) {
          EXPECT_EQ(cells.mIntA(i, d), static_cast<Int>(j + d));
      EXPECT_NUM_EQ(cells.mNumA(i, d), static_cast<Num>(j) / 2 + d);
    }
    ++i;
  }
}

TEST(fixed_container_test, iterator_sort) {
  FC2D cells(100);
  cells.push_cell(10);
  init_variables(cells);
  check_variables(cells, 0, cells.size());

  auto check_inversely_sorted = [&]() {
    for (Ind i = 0, e = 10; i < e; ++i) {
          EXPECT_EQ(cells.mInt(i), static_cast<Int>(9 - i));
      EXPECT_NUM_EQ(cells.mNum(i), static_cast<Num>(9 - i) / 2);
      for (SInd d = 0; d < FC2D::nd; ++d) {
            EXPECT_EQ(cells.mIntA(i, d), static_cast<Int>(9 - i + d));
        EXPECT_NUM_EQ(cells.mNumA(i, d), static_cast<Num>(9 - i) / 2 + d);
      }
    }
  };

  algorithm::sort(cells, [](FC2DRef i, FC2DRef j) {
    return i.mNum() > j.mNum();
  });

  check_inversely_sorted();
}


/// Test Cells with variable number of nodes:
template<class C> void plotCellNodes2D(C& cells) {
  TRACE_IN_();
  for (Ind cIdx = 0, endCIdx = cells.size(); cIdx < endCIdx; ++cIdx) {
    DBG_ON("cId:", cIdx, "[nB,nE]: [", cells.first_node(cIdx),        \
        ",", cells.last_node(cIdx), "]", "#nodes:", cells.node_size(cIdx));
    plotCellNodeVariables2D(cells, cIdx);
  }
  TRACE_OUT();
}

template<class C> void plotCellNodeVariables2D(C& cells, Ind cIdx) {
  TRACE_IN((cIdx));
  std::cerr << "# cId: " << cIdx << " ";
  for (Ind nIdx = cells.first_node(cIdx), endNodeId = cells.last_node(cIdx);
       nIdx < endNodeId; ++nIdx) {
    std::cerr << " | i:" << cells.mIntN(nIdx) <<  " n:" << cells.mNumNE(nIdx);
  }
  std::cerr << std::endl;
  TRACE_OUT();
}

template<class C> void initNodalVariables2D(C& cells) {
  TRACE_IN_();
  init_variables(cells);

  for (Ind cIdx = 0, endCIdx = cells.size(); cIdx < endCIdx; ++cIdx) {
    for (Ind nIdx = cells.first_node(cIdx),
         endNodeId = cells.last_node(cIdx); nIdx < endNodeId; ++nIdx) {
      cells.mIntN(nIdx) = nIdx;
      cells.mNumN(nIdx) = 2*nIdx;
    }
  }
  TRACE_OUT();
}

template<class C> void checkNodalVariables2D
(C& cells, C& cells_copy, Ind fromCIdx, Ind toCIdx) {
  TRACE_IN_();
  check_variables(cells, fromCIdx, toCIdx);

  for (Ind cIdx = fromCIdx, endCIdx = toCIdx; cIdx < endCIdx; ++cIdx) {
    for (Ind ni = cells.first_node(cIdx), nei = cells.last_node(cIdx),
            nc = cells_copy.first_node(cIdx), nec = cells_copy.last_node(cIdx);
            nc < nec; ++nc, ++ni) {
      EXPECT_EQ((nei - ni), (nec - nc));
      EXPECT_EQ(cells.mIntN(ni), cells_copy.mIntN(nc));
      EXPECT_NUM_EQ(cells.mNumN(ni), cells_copy.mNumN(nc));
    }
  }
  TRACE_OUT();
}

template<class C> void copy_last_half_to_beginning(C& cells, C& cells_copy) {
  TRACE_IN_();
  /// Copy cells:
  for (Ind i = cells.size() / 2, e = cells.size(), j = cells.size() / 2;
       i < e; ++i) {
    cells.copy_cell(i, i - j);
  }
  /// Test:
  /// First half:
  for (Ind cIdx = 0, endCIdx = cells.size() / 2; cIdx < endCIdx; ++cIdx) {
    EXPECT_EQ(cells.mInt(cIdx), cells_copy.mInt(cIdx + endCIdx));
    EXPECT_NUM_EQ(cells.mNum(cIdx), cells_copy.mNum(cIdx + endCIdx));
    for (SInd d = 0; d < C::nd; ++d) {
      EXPECT_EQ(cells.mIntA(cIdx, d), cells_copy.mIntA(cIdx + endCIdx, d));
      EXPECT_NUM_EQ(cells.mNumA(cIdx, d), cells_copy.mNumA(cIdx + endCIdx, d));
    }
    for (Ind ni = cells.first_node(cIdx), nei = cells.last_node(cIdx),
             nc = cells_copy.first_node(cIdx + endCIdx),
            nec = cells_copy.last_node(cIdx + endCIdx);
        nc < nec; ++nc, ++ni) {
      EXPECT_TRUE(ni < nei);
      EXPECT_EQ(cells.mIntN(ni), cells_copy.mIntN(nc));
      EXPECT_NUM_EQ(cells.mNumN(ni), cells_copy.mNumN(nc));
    }
  }

  /// Second half.
  checkNodalVariables2D(cells, cells_copy, cells.size() / 2, cells.size());
  TRACE_OUT();
}

NC2D create_linear_increasing_nodal_2d_container(Ind noCells, Ind maxNoCells) {
  Ind linearMaxNoNodes = noCells *(noCells + 1) / 2;

  NC2D cells(maxNoCells, 2 * linearMaxNoNodes);
  for (Ind cIdx = 0; cIdx < noCells; ++cIdx)
    cells.push_cell(1, cIdx + 1);
  initNodalVariables2D(cells);

  EXPECT_EQ(cells.size(), noCells);
  EXPECT_EQ(cells.node_size(), linearMaxNoNodes);
  return cells;
}

NC2D create_linear_decreasing_nodal_2d_container(Ind noCells, Ind maxNoCells) {
  Ind linearMaxNoNodes = noCells * (noCells + 1) / 2;

  NC2D cells(maxNoCells, 2 * linearMaxNoNodes);
  for (Ind cIdx = noCells; cIdx > 0; --cIdx) {
    cells.push_cell(1, cIdx);
  }
  initNodalVariables2D(cells);

  EXPECT_EQ(cells.size(), noCells);
  EXPECT_EQ(cells.node_size(), linearMaxNoNodes);
  return cells;
}

template<class C> Ind erase_even_cells_and_verify(C& cells, C& cells_copy) {
  TRACE_IN_();
  /// Erase even cells:
  auto oldSize = cells.size();
  auto pred = [](Ind cIdx) { return cIdx % 2 == 0; };
  cells.erase_remove_if(pred);
  /// Check:
  for (NC2D::cell_size_type j = 0, ej = oldSize, i = 0; j < ej; ++j) {
    if (j % 2 == 0) continue;
    /// cell values:
    EXPECT_EQ(cells.mInt(i), cells_copy.mInt(j));
    EXPECT_NUM_EQ(cells.mNum(i), cells_copy.mNum(j));
    for (SInd d = 0; d < C::nd; ++d) {
      EXPECT_EQ(cells.mIntA(i, d), cells_copy.mIntA(j, d));
      EXPECT_NUM_EQ(cells.mNumA(i, d), cells_copy.mNumA(j, d));
    }
    /// nodal values:
    for (Ind ni = cells.first_node(i), nei = cells.last_node(i),
             nc = cells_copy.first_node(j), nec = cells_copy.last_node(j);
         nc < nec; ++nc, ++ni) {
      EXPECT_TRUE(ni < nei);
      EXPECT_EQ(cells.mIntN(ni), cells_copy.mIntN(nc));
      EXPECT_NUM_EQ(cells.mNumN(ni), cells_copy.mNumN(nc));
    }
    ++i;
  }
  TRACE_OUT();
  return cells.size();
}

/// Nodal container where all cells have the same #of nodes:
TEST(variable_container_test, constant_copy) {
  Ind noCells = 10, maxNoCells = 10, maxNoNodes = 100;

  /// Constant container:
  Ind nodesPerCell = 10;

  NC2D cells(maxNoCells, maxNoNodes);
  cells.push_cell(noCells, nodesPerCell);
  initNodalVariables2D(cells);

  NC2D cells_copy(cells);

  checkNodalVariables2D(cells, cells_copy, 0, cells.size());

  copy_last_half_to_beginning(cells, cells_copy);
}

TEST(variable_container_test, constant_erase) {
  Ind noCells = 10, maxNoCells = 10, maxNoNodes = 100;

  /// Constant container:
  Ind nodesPerCell = 10;

  /// Create container and copy:
  NC2D cells(maxNoCells, maxNoNodes);
  cells.push_cell(noCells, nodesPerCell);
  initNodalVariables2D(cells);
  NC2D cells_copy(cells);
  checkNodalVariables2D(cells, cells_copy, 0, cells.size());

  /// Erase even cells:
  auto size = erase_even_cells_and_verify(cells, cells_copy);

  /// Check size:
  EXPECT_EQ(size, maxNoCells / 2);
}

/// Nodal container where the #of nodes in the cells is incremented linearly
TEST(variable_container_test, linear_incr_copy) {
  Ind noCells = 10, maxNoCells = 10;  // maxNoNodes = 100;

  NC2D cells(create_linear_increasing_nodal_2d_container(noCells, maxNoCells));
  NC2D cells_copy(cells);
  checkNodalVariables2D(cells, cells_copy, 0, cells.size());

  copy_last_half_to_beginning(cells, cells_copy);
}

TEST(variable_container_test, linear_incr_erase) {
  Ind noCells = 10, maxNoCells = 10;  // maxNoNodes = 100;

  NC2D cells(create_linear_increasing_nodal_2d_container(noCells, maxNoCells));
  NC2D cells_copy(cells);
  checkNodalVariables2D(cells, cells_copy, 0, cells.size());

  auto size = erase_even_cells_and_verify(cells, cells_copy);

  /// Check size:
  EXPECT_EQ(size, maxNoCells / 2);
  EXPECT_EQ(cells.size(), static_cast<Ind>(maxNoCells / 2));
}

/// Nodal container where the #of nodes in the cells is decremented linearly
TEST(variable_container_test, linear_decr_copy) {
  Ind noCells = 10, maxNoCells = 10;  // maxNoNodes = 100;

  NC2D cells(create_linear_decreasing_nodal_2d_container(noCells, maxNoCells));
  NC2D cells_copy(cells);
  checkNodalVariables2D(cells, cells_copy, 0, cells.size());

  copy_last_half_to_beginning(cells, cells_copy);
}

TEST(variable_container_test, linear_decr_erase) {
  Ind noCells = 10, maxNoCells = 10;  // maxNoNodes = 100;

  NC2D cells(create_linear_decreasing_nodal_2d_container(noCells, maxNoCells));
  NC2D cells_copy(cells);
  checkNodalVariables2D(cells, cells_copy, 0, cells.size());

  auto size = erase_even_cells_and_verify(cells, cells_copy);

  /// Check size:
  EXPECT_EQ(size, maxNoCells / 2);
  EXPECT_EQ(cells.size(), static_cast<Ind>(maxNoCells / 2));
}


TEST(variable_container_test, single_index_storage) {
  Ind maxNoCells = 5, maxNoNodes = 20;

  NC2D cells(maxNoCells, maxNoNodes);
  using nst = NC2D::node_size_type;
  cells.push_cell(1, 1);  cells.push_cell(1, 2);  cells.push_cell(1, 3);

  // plotCellNodes2D(cells);
  DBG("");
  EXPECT_EQ(cells.node_size(0), nst{1});
  EXPECT_EQ(cells.node_size(1), nst{2});
  EXPECT_EQ(cells.node_size(2), nst{3});

  DBG("cells.erase_remove_if([&cells](Ind i) { return i == 2; });");
  cells.erase_remove_if([&cells](Ind i) { return i == 2; });
  // plotCellNodes2D(cells);
  EXPECT_EQ(cells.node_size(0), nst{1});
  EXPECT_EQ(cells.node_size(1), nst{2});
  DBG("");

  DBG("cells.push_back(1, 5);");
  cells.push_cell(1, 5);
  // plotCellNodes2D(cells);
  EXPECT_EQ(cells.node_size(0), nst{1});
  EXPECT_EQ(cells.node_size(1), nst{2});
  EXPECT_EQ(cells.node_size(2), nst{5});
  DBG("");

  DBG("cells.erase_remove_if([&cells](Ind i) { return i == 0; });");
  cells.erase_remove_if([&cells](Ind i) { return i == 0; });
  // plotCellNodes2D(cells);
  DBG("");

  DBG("cells.push_back(1, 3);");
  cells.push_cell(1, 3);
  // plotCellNodes2D(cells);
  DBG("");

  DBG("cells.erase_remove_if([&cells](Ind i) { return i == 2; });");
  cells.erase_remove_if([&cells](Ind i) { return i == 2; });
  // plotCellNodes2D(cells);
  DBG("");

  DBG("cells.push_back(1, 4);");
  cells.push_cell(1, 4);
  // plotCellNodes2D(cells);
  DBG("");
}

////////////////////////////////////////////////////////////////////////////////
#undef ENABLE_DBG_
////////////////////////////////////////////////////////////////////////////////
