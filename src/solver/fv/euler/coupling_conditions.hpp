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

/// \brief Euler-equations coupling conditions
namespace coupling {

/// \brief Coupling condition with another Euler solver
template<class EulerSolver>
struct Euler : fv::bc::Condition<Euler<EulerSolver>> {
  Euler(EulerSolver& thisSolver, const EulerSolver& otherSolver)
      : s(thisSolver), os_(otherSolver) {}

  /// \brief Applies the boundary condition to a range of ghost cells
  template<class _>
  void operator()(_, const CellIdx thisGhostIdx) const noexcept {
      CellIdx thisBndryIdx, otherBndryIdx;
      std::tie(thisBndryIdx, otherBndryIdx)
        = fv::coupling::local_bndry_ids(thisGhostIdx, s, os_);
      s.Q(_(), thisGhostIdx) = os_.Q(_(), otherBndryIdx);
  }

  EulerSolver& s;
  EulerSolver const& os_;
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
