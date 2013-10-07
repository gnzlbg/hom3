#ifndef HOM3_SOLVERS_FV_HEAT_COUPLING_CONDITIONS_HPP_
#define HOM3_SOLVERS_FV_HEAT_COUPLING_CONDITIONS_HPP_
////////////////////////////////////////////////////////////////////////////////
/// \file \brief Implements the Navier-Stokes coupling conditions.
////////////////////////////////////////////////////////////////////////////////
#include "solver/fv/heat/boundary_conditions.hpp"
#include "solver/fv/coupling.hpp"
////////////////////////////////////////////////////////////////////////////////
namespace hom3 { namespace solver { namespace fv { namespace heat {
////////////////////////////////////////////////////////////////////////////////
namespace bc {
////////////////////////////////////////////////////////////////////////////////

/// \brief Coupling conditions
namespace coupling {

/// \brief Coupling condition with a CNS solver
template<class HeatSolver, class CNSSolver> struct CNS : Dirichlet<HeatSolver> {
  CNS(HeatSolver& heatSolver, const CNSSolver& cnsSolver)
    : Dirichlet<HeatSolver>(heatSolver, [&](const CellIdx ghostIdx) {
        return surface_temperature(ghostIdx, heatSolver, cnsSolver);})
  {}

  inline Num surface_temperature
  (const CellIdx heatGhostIdx, const HeatSolver& hs_,
    const CNSSolver& cns_) const noexcept {
    CellIdx heatBndryIdx, cnsBndryIdx;
    std::tie(heatBndryIdx, cnsBndryIdx)
        = fv::coupling::local_bndry_ids(heatGhostIdx, hs_, cns_);

    const auto T_solid = hs_.T(lhs, heatBndryIdx);  // T from last time step
    const auto T_fluid = cns_.T(lhs, cnsBndryIdx);  // T from last time step

    const auto T_srfc = 0.5 * (T_fluid + T_solid);  // \todo interpolate better!

    return T_srfc;
  }
};

/// \brief Coupling condition with another Heat-equation solver
template<class HeatSolver> struct Heat : Dirichlet<HeatSolver> {
  Heat(HeatSolver& thisSolver, const HeatSolver& otherSolver)
    : Dirichlet<HeatSolver>(thisSolver, [&](const CellIdx ghostIdx) {
        return surface_temperature(ghostIdx, thisSolver, otherSolver);})
  {}

  inline Num surface_temperature
  (const CellIdx thisGhostIdx, const HeatSolver& ts_,
    const HeatSolver& os_) const noexcept {
    CellIdx thisBndryIdx, otherBndryIdx;
    std::tie(thisBndryIdx, otherBndryIdx)
        = fv::coupling::local_bndry_ids(thisGhostIdx, ts_, os_);

    const auto T_this = ts_.T(lhs, thisBndryIdx);  // T from current RK step
    const auto T_other = os_.T(lhs, otherBndryIdx);  // T from last time step

    const auto T_srfc = 0.5 * (T_this + T_other);  // \todo interpolate better!

    return T_srfc;
  }
};


}  // namespace coupling

////////////////////////////////////////////////////////////////////////////////
}  // namespace bc
}  // namespace heat
}  // namespace fv
}  // namespace solver
}  // namespace hom3
////////////////////////////////////////////////////////////////////////////////
#endif
