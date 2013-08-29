#ifndef HOM3_SOLVERS_FV_EULER_COUPLING_CONDITIONS_HPP_
#define HOM3_SOLVERS_FV_EULER_COUPLING_CONDITIONS_HPP_
////////////////////////////////////////////////////////////////////////////////
/// \file \brief Implements the Euler-equations coupling conditions.
////////////////////////////////////////////////////////////////////////////////
#include "solver/fv/coupling.hpp"
////////////////////////////////////////////////////////////////////////////////
namespace hom3 { namespace solver { namespace fv { namespace euler {
////////////////////////////////////////////////////////////////////////////////
namespace bc {
////////////////////////////////////////////////////////////////////////////////

/// \brief Coupling conditions
namespace coupling {

/// \brief Coupling condition with another Euler solver Heat solver
template<class EulerSolver>
struct Euler : Dirichlet<EulerSolver> {
  Euler(EulerSolver& thisSolver, const EulerSolver& otherSolver)
      : Dirichlet<EulerSolver>(thisSolver, [&](const CellIdx bndryIdx,
                                               const SInd v) {
          return surface_variables(bndryIdx, v, thisSolver, otherSolver);})
  {}

  inline Num surface_variables
  (const CellIdx thisGhostIdx, const SInd v, const EulerSolver& ts_,
   const EulerSolver& os_) const noexcept {
    CellIdx thisBndryIdx, otherBndryIdx;
    std::tie(thisBndryIdx, otherBndryIdx)
      = fv::coupling::local_bndry_ids(thisGhostIdx, ts_, os_);

    return os_.Q(lhs, otherBndryIdx, v);
  }
};

}  // namespace coupling

////////////////////////////////////////////////////////////////////////////////
}  // namespace bc
}  // namespace euler
}  // namespace fv
}  // namespace solver
}  // namespace hom3
////////////////////////////////////////////////////////////////////////////////
#endif
