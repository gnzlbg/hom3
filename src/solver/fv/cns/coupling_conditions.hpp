#ifndef HOM3_SOLVERS_FV_CNS_COUPLING_CONDITIONS_HPP_
#define HOM3_SOLVERS_FV_CNS_COUPLING_CONDITIONS_HPP_
////////////////////////////////////////////////////////////////////////////////
/// \file \brief Implements the Navier-Stokes coupling conditions.
////////////////////////////////////////////////////////////////////////////////
#include "solver/fv/cns/boundary_conditions.hpp"
#include "solver/fv/cns/moving_boundary_conditions.hpp"
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

template<class CNSSolver, class HeatSolver>
struct MBHeat : moving_wall::IsothermalNoSlip<CNSSolver> {
  template<class G, class V, class A>
  MBHeat(CNSSolver& cnsSolver, const HeatSolver& heatSolver, G&& g, V&& v, A&& a)
      : moving_wall::IsothermalNoSlip<CNSSolver>(cnsSolver, g, v, a,
                                                [&](const CellIdx ghostIdx) {
        return surface_temperature(ghostIdx, cnsSolver, heatSolver);})
  {}

  inline Num surface_temperature
  (const CellIdx cnsGhostIdx, const CNSSolver& cns_,
    const HeatSolver& hs_) const noexcept {

    CellIdx cnsBndryIdx = cnsGhostIdx;
    CellIdx heatBndryIdx = cns_.grid().cell_idx
                          (cns_.node_idx(cnsBndryIdx), hs_.solver_idx());
    ASSERT(is_valid(heatBndryIdx), "!!");
    ASSERT(is_valid(cnsBndryIdx), "!!");

    /// find closest cell in heat domain
    auto in = hs_.find_interpolation_neighbors(heatBndryIdx);
    auto nghbrs = memory::stack::vector<CellIdx, 82>{};
    nghbrs.reserve(40);
    Num dist = std::numeric_limits<Num>::max();
    CellIdx heatInternalIdx = invalid<CellIdx>();
    Int ns = 0;
    for (auto n : in) {
      if (is_valid(n) && hs_.exists(n) && is_valid(hs_.node_idx(n))) {
        nghbrs.push_back(n);
        auto tmp =  hs_.find_interpolation_neighbors(n);
        for (auto t : tmp) {
          if (is_valid(t) && hs_.exists(t) && is_valid(hs_.node_idx(t))) {
            nghbrs.push_back(t);
          }
        }
      }
    }

    for (auto n : nghbrs) {
      if (is_valid(n) && !hs_.is_ghost_cell(n)
          && !is_valid(hs_.cut_by_which_moving_boundary(n))) {
        ++ns;
        if (hs_.cell_dx(heatBndryIdx, n) < dist) {
          dist = hs_.cell_dx(heatBndryIdx, n);
          heatInternalIdx = n;
        }
      }
    }

    if(!is_valid(heatInternalIdx)) {
      std::cerr << "hIdx: " << heatBndryIdx << " cnsIdx: " << cnsBndryIdx
                << " no_in: " << in.size()
                << " no_valid_in: " << nghbrs.size() << std::endl;
      std::cerr << "ns: ";
      for (auto n : nghbrs) {
        std::cerr << n << " ";
      }
      std::cerr << std::endl;
      TERMINATE("HEAT CC: INVALID IDX!");
    }
    Num T_srfc = hs_.T(lhs, heatInternalIdx);
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
