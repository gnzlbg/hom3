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

/// \brief Coupling condition with another Euler solver
template<class EulerSolver>
struct Euler : fv::bc::Condition<Euler<EulerSolver>> {
  Euler(EulerSolver& thisSolver, const EulerSolver& otherSolver)
      : s(thisSolver), os_(otherSolver) {}

  /// \brief Applies the boundary condition to a range of ghost cells
  template<class _>
  void operator()(_, Range<CellIdx> ghost_cells) const noexcept {
    for (auto thisGhostIdx : ghost_cells) {
      CellIdx thisBndryIdx, otherBndryIdx;
      std::tie(thisBndryIdx, otherBndryIdx)
        = fv::coupling::local_bndry_ids(thisGhostIdx, s, os_);
      for (SInd v = 0; v < s.nvars; ++v) {
        s.Q(_(), thisGhostIdx, v) = os_.Q(_(), otherBndryIdx, v);
      }
    }
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
