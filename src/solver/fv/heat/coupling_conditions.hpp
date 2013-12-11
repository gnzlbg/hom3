#ifndef HOM3_SOLVERS_FV_HEAT_COUPLING_CONDITIONS_HPP_
#define HOM3_SOLVERS_FV_HEAT_COUPLING_CONDITIONS_HPP_
////////////////////////////////////////////////////////////////////////////////
/// \file \brief Implements the Navier-Stokes coupling conditions.
////////////////////////////////////////////////////////////////////////////////
#include "solver/fv/heat/boundary_conditions.hpp"
#include "solver/fv/heat/moving_boundary_conditions.hpp"
#include "solver/fv/coupling.hpp"
////////////////////////////////////////////////////////////////////////////////
namespace hom3 { namespace solver { namespace fv { namespace heat {
////////////////////////////////////////////////////////////////////////////////
namespace bc {
////////////////////////////////////////////////////////////////////////////////

/// \brief Heat equation coupling conditions
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

/// \brief Coupling condition with a CNS solver
template<class HeatSolver, class CNSSolver> struct MBCNS : mb::Dirichlet<HeatSolver> {
template<class G, class V>
MBCNS(HeatSolver& heatSolver, const CNSSolver& cnsSolver, G&& geometry, V&& vel)
    : mb::Dirichlet<HeatSolver>(heatSolver, geometry, vel,
                                [&](const CellIdx ghostIdx) {
    return surface_temperature(ghostIdx, heatSolver, cnsSolver);})
  {}

  inline Num surface_temperature
  (const CellIdx heatGhostIdx, const HeatSolver& hs_,
    const CNSSolver& cns_) const noexcept {
    CellIdx heatBndryIdx = heatGhostIdx;
    CellIdx cnsBndryIdx = hs_.grid().cell_idx
                          (hs_.node_idx(heatBndryIdx), cns_.solver_idx());
    ASSERT(is_valid(heatBndryIdx), "!!");
    ASSERT(is_valid(cnsBndryIdx), "!!");

    /// find closest cell in cns domain
    auto in = cns_.find_interpolation_neighbors(cnsBndryIdx);
    auto nghbrs = memory::stack::vector<CellIdx, 82>{};
    nghbrs.reserve(40);
    Num dist = std::numeric_limits<Num>::max();
    CellIdx cnsInternalIdx = invalid<CellIdx>();
    Int ns = 0;
    for (auto n : in) {
      if (is_valid(n) && cns_.exists(n) && is_valid(cns_.node_idx(n))) {
        nghbrs.push_back(n);
        auto tmp =  cns_.find_interpolation_neighbors(n);
        for (auto t : tmp) {
          if (is_valid(t) && cns_.exists(t) && is_valid(cns_.node_idx(t))) {
            nghbrs.push_back(t);
          }
        }
      }
    }

    for (auto n : nghbrs) {
      if (is_valid(n) && !cns_.is_ghost_cell(n)
          && !is_valid(cns_.cut_by_which_moving_boundary(n))) {
        ++ns;
        if (cns_.cell_dx(cnsBndryIdx, n) < dist) {
          dist = cns_.cell_dx(cnsBndryIdx, n);
          cnsInternalIdx = n;
        }
      }
    }

    if(!is_valid(cnsInternalIdx)) {
      std::cerr << "hIdx: " << heatGhostIdx << " cnsIdx: " << cnsBndryIdx
                << " no_in: " << in.size()
                << " no_valid_in: " << nghbrs.size() << std::endl;
      std::cerr << "ns: ";
      for (auto n : nghbrs) {
        std::cerr << n << " ";
      }
      std::cerr << std::endl;
      TERMINATE("HEAT CC: INVALID IDX!");
    }

    Num T_srfc = cns_.T(lhs, cnsInternalIdx);

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
