#ifndef HOM3_CONTAINERS_SEQUENTIAL_TESTS_HELPERS_HPP_
#define HOM3_CONTAINERS_SEQUENTIAL_TESTS_HELPERS_HPP_
////////////////////////////////////////////////////////////////////////////////
/// \file \brief Helper functions for containers::sequential unit-tests
////////////////////////////////////////////////////////////////////////////////
/// Includes:
#include "containers/sequential/tests/fixedcell.hpp"
#include "containers/sequential/tests/nodalcell.hpp"
/// Options:
#define ENABLE_DBG_ 0
#include "misc/dbg.hpp"
////////////////////////////////////////////////////////////////////////////////
namespace hom3 { namespace container { namespace sequential {
////////////////////////////////////////////////////////////////////////////////

namespace test {

using FC2D    = fixed::Container<2>;
using FC3D    = fixed::Container<3>;
using FC2DRef = FC2D::reference;
using FC3DRef = FC3D::reference;
using FC2DVal = FC2D::value_type;
using FC3DVal = FC3D::value_type;

using NC2D    = nodal::Container<2>;
using NC3D    = nodal::Container<3>;

template<SInd nd>
void init_cell(typename fixed::Container<nd>::reference&& c, Ind incr) {
  c.mInt() = 0 + static_cast<Int>(incr);
  c.mNum() = 1 + static_cast<Num>(incr);
  for (SInd d = 0; d < nd; ++d) {
    c.mIntA(d) = d + static_cast<Int>(incr);
    c.mNumA(d) = d + static_cast<Num>(incr);
  }
}

template<SInd nd>
void check_cell(typename fixed::Container<nd>::value_type&& c, Ind incr) {
      EXPECT_EQ(c.mInt(), 0 + static_cast<Int>(incr));
  EXPECT_NUM_EQ(c.mNum(), 1 + static_cast<Num>(incr));
  for (SInd d = 0; d < nd; ++d) {
        EXPECT_EQ(c.mIntA(d), d + static_cast<Int>(incr));
    EXPECT_NUM_EQ(c.mNumA(d), d + static_cast<Num>(incr));
  }
}

template<SInd nd>
void check_cell(typename fixed::Container<nd>::reference&& c, Ind incr) {
      EXPECT_EQ(c.mInt(), 0 + static_cast<Int>(incr));
  EXPECT_NUM_EQ(c.mNum(), 1 + static_cast<Num>(incr));
  for (SInd d = 0; d < nd; ++d) {
        EXPECT_EQ(c.mIntA(d), d + static_cast<Int>(incr));
    EXPECT_NUM_EQ(c.mNumA(d), d + static_cast<Num>(incr));
  }
}

template<SInd nd, class T1, class T2>
void check_cell_equality(T1&& c1, T2&& c2) {
      EXPECT_EQ(c1.mInt(), c2.mInt());
  EXPECT_NUM_EQ(c1.mNum(), c2.mNum());
  for (SInd d = 0; d < nd; ++d) {
        EXPECT_EQ(c1.mIntA(d), c1.mIntA(d));
    EXPECT_NUM_EQ(c1.mNumA(d), c2.mNumA(d));
  }
}

template<class C> void init_variables(C& cells) {
  TRACE_IN_();
  for (Ind i = 0, e = cells.size(); i < e; ++i) {
    cells.mInt(i) = static_cast<Int>(i);
    cells.mNum(i) = static_cast<Num>(i) / 2;
    for (SInd d = 0; d < C::nd; ++d) {
      cells.mIntA(i, d) = static_cast<Int>(i + d);
      cells.mNumA(i, d) = static_cast<Num>(i) / 2 + d;
    }
  }
  TRACE_OUT();
}

template<class C> void check_variables(C& cells, Ind fromId, Ind toId) {
  TRACE_IN_();
  for (Ind i = fromId, e = toId; i < e; ++i) {
    EXPECT_EQ(cells.mInt(i)    , static_cast<Int>(i));
    EXPECT_NUM_EQ(cells.mNum(i), static_cast<Num>(i) / 2);
    for (SInd d = 0; d < C::nd; ++d) {
      EXPECT_EQ(cells.mIntA(i, d)    , static_cast<Int>(i + d));
      EXPECT_NUM_EQ(cells.mNumA(i, d), static_cast<Num>(i) / 2 + d);
    }
  }
  TRACE_OUT();
}

template<class C> void check_variables(C& cells) {
  check_variables(cells, cells.first(), cells.size());
}

template<class CellRef>
void check_vars(CellRef cell) {
  EXPECT_EQ(cell.mInt()    , static_cast<Int>(cell.index()));
  EXPECT_NUM_EQ(cell.mNum(), static_cast<Num>(cell.index()) / 2);
  for (SInd d = 0; d < CellRef::nd; ++d) {
    EXPECT_EQ(cell.mIntA(d)    , static_cast<Int>(cell.index() + d));
    EXPECT_NUM_EQ(cell.mNumA(d), static_cast<Num>(cell.index()) / 2 + d);
  }
};

}  // namespace test

////////////////////////////////////////////////////////////////////////////////
}  // namespace sequential
}  // namespace container
}  // namespace hom3
////////////////////////////////////////////////////////////////////////////////
#undef ENABLE_DBG_
////////////////////////////////////////////////////////////////////////////////
#endif
