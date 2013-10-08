#ifndef HOM3_SOLVERS_FV_CNS_COUPLING_CONDITIONS_HPP_
#define HOM3_SOLVERS_FV_CNS_COUPLING_CONDITIONS_HPP_
////////////////////////////////////////////////////////////////////////////////
/// \file \brief Implements the Navier-Stokes coupling conditions.
////////////////////////////////////////////////////////////////////////////////
#include "solver/fv/cns/boundary_conditions.hpp"
#include "solver/fv/coupling.hpp"
////////////////////////////////////////////////////////////////////////////////
namespace hom3 { namespace solver { namespace fv { namespace cns {
////////////////////////////////////////////////////////////////////////////////
namespace bc {
////////////////////////////////////////////////////////////////////////////////

/// \brief Navier-Stokes coupling conditions
namespace coupling {

/// \brief Coupling condition with a Heat solver
///
/// Isothermal no-slip wall with surface temperature
/// $T_\mathrm{srfc} = \frac{T_\mathrm{f} + T_\mathrm{s}}{2}$
template<class CNSSolver, class HeatSolver>
struct Heat : wall::IsothermalNoSlip<CNSSolver> {
  Heat(CNSSolver& cnsSolver, const HeatSolver& heatSolver)
    : wall::IsothermalNoSlip<CNSSolver>(cnsSolver, [&](const CellIdx ghostIdx) {
        return surface_temperature(ghostIdx, cnsSolver, heatSolver);})
  {}

  inline Num surface_temperature
  (const CellIdx cnsGhostIdx, const CNSSolver& cns_,
    const HeatSolver& hs_) const noexcept {
    CellIdx cnsBndryIdx, heatBndryIdx;
    std::tie(cnsBndryIdx, heatBndryIdx)
      = fv::coupling::local_bndry_ids(cnsGhostIdx, cns_, hs_);

    // const auto T_fluid = cns_.T(rhs, cnsBndryIdx);  // T from current RK step
    const auto T_fluid = cns_.T(lhs, cnsBndryIdx);  // T from current RK step
    const auto T_solid = hs_.T(lhs, heatBndryIdx);  // T from last time step

    const auto T_srfc = 0.5 * (T_fluid + T_solid);  // \todo interpolate better!

    return T_srfc;
  }
};

}  // namespace coupling

////////////////////////////////////////////////////////////////////////////////
}  // namespace bc
}  // namespace cns
}  // namespace fv
}  // namespace solver
}  // namespace hom3
////////////////////////////////////////////////////////////////////////////////
#endif
