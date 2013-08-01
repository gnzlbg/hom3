#ifndef HOM3_MISC_INDICES_HPP_
#define HOM3_MISC_INDICES_HPP_
////////////////////////////////////////////////////////////////////////////////
#include "types.hpp"
#include "integer.hpp"
#include "ranges.hpp"
////////////////////////////////////////////////////////////////////////////////

/// \name Indexing Types
///@{

namespace detail_ {
struct node_idx_tag {};
struct cell_idx_tag {};
struct solver_idx_tag {};
}

using NodeIdx   = Integer<Ind,detail_::node_idx_tag>;
using CellIdx   = Integer<Ind,detail_::cell_idx_tag>;
using SolverIdx = Integer<SInd,detail_::solver_idx_tag>;
template<SInd nRows>  using NodeIdxA  = EigenStaticVector<NodeIdx,nRows>;
template<SInd nd = 1> using NodeIdxM  = EigenColMajor<NodeIdx, nd>;
template<SInd nd = 1> using NodeIdxRM = EigenRowMajor<NodeIdx, nd>;
template<SInd nRows>  using CellIdxA  = EigenStaticVector<CellIdx,nRows>;
template<SInd nd = 1> using CellIdxM  = EigenColMajor<CellIdx, nd>;
template<SInd nd = 1> using CellIdxRM = EigenRowMajor<CellIdx, nd>;

/// \brief Cast range of to to range of primitive types
template<class T>
inline auto primitive_cast() -> RangeTransformer<T, typename T::value_type> {
 return {[&](const T i){ return i(); }};
}

///@}

////////////////////////////////////////////////////////////////////////////////
#endif
