#ifndef HOM3_GEOMETRY_CELL_CARTESIAN_HPP_
#define HOM3_GEOMETRY_CELL_CARTESIAN_HPP_
////////////////////////////////////////////////////////////////////////////////
#include "globals.hpp"
//////////////////////////////////////////////////////////////////////////

namespace hom3 { namespace geometry { namespace cell { namespace cartesian {

template<SInd nd> constexpr Num volume(const Num dx) {
  return dx * volume<nd-1>(dx);
}
template<> constexpr Num volume<1>(const Num dx) { return dx; }


}}}} // hom3::geometry::cell::cartesian namespace

////////////////////////////////////////////////////////////////////////////////
#endif
