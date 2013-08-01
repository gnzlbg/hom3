#ifndef HOM3_SOLVERS_FV_COUPLING_HPP_
#define HOM3_SOLVERS_FV_COUPLING_HPP_
////////////////////////////////////////////////////////////////////////////////
/// Options:
#define ENABLE_DBG_ 0
#include "../../misc/dbg.hpp"
////////////////////////////////////////////////////////////////////////////////

namespace solver { namespace fv { namespace coupling {

/// \brief Computes the local boundary cell ids of solver0 and solver1
///
/// \returns (solver0LocalBdryCellId,solver1LocalBndryCellId)
template<class Solver0, class Solver1>
std::tuple<CellIdx,CellIdx> local_bndry_ids(const CellIdx solver0GhostCellId,
                                            Solver0& solver0, const Solver1& solver1) {

  CellIdx solver0BndryCellId; SInd solver0BndryCellPos;
  std::tie(solver0BndryCellId,solver0BndryCellPos)
      = solver0.boundary_cell_id(solver0GhostCellId);

  auto solver0GhostCellPos
      = solver0.grid().nodes().opposite_nghbr_position(solver0BndryCellPos);

  auto solver0GlobalBndryCellId = solver0.node_idx(solver0BndryCellId);
  auto solver1GlobalBndryCellId
      = solver0.grid().nodes().find_samelvl_nghbr(solver0GlobalBndryCellId,
                                                        solver0GhostCellPos);
  auto solver1LocalBndryCellId
      = solver0.grid().nodes().cell_id(solver1GlobalBndryCellId,
                                             solver1.solver_idx());
  return std::make_tuple(solver0BndryCellId, solver1LocalBndryCellId);
}

namespace bc {

template<class HeatSolver, class EulerSolver>
struct HeatEulerWall : boundary::Condition {
  HeatEulerWall(HeatSolver& heatSolver, const EulerSolver& eulerSolver)
      : heatSolver_(heatSolver), eulerSolver_(eulerSolver) { }

  void apply(AnyRange<CellIdx> heatGhostCells) const override {
    for(auto heatGhostCellId : heatGhostCells) {
      CellIdx heatBndryCellId, eulerBndryCellId;
      std::tie(heatBndryCellId, eulerBndryCellId)
          = local_bndry_ids(heatGhostCellId,heatSolver_,eulerSolver_);

      const Num surfaceEnergy = 0.5 * (heatSolver_.cells().lhs(heatBndryCellId,0)
                                 + eulerSolver_.cells().lhs(eulerBndryCellId,
                                                             EulerSolver::V::E()));
      heatSolver_.cells().lhs(heatGhostCellId,0)
          = 2.0 * surfaceEnergy - heatSolver_.cells().lhs(heatBndryCellId,0);
      heatSolver_.cells().rhs(heatGhostCellId,0)
          = 2.0 * surfaceEnergy - heatSolver_.cells().rhs(heatBndryCellId,0);
    }
  }

 private:
  HeatSolver& heatSolver_;
  const EulerSolver& eulerSolver_;
};


template<class EulerSolver, class HeatSolver>
struct EulerHeatWall : boundary::Condition {

  EulerHeatWall(EulerSolver& eulerSolver, const HeatSolver& heatSolver)
      : eulerSolver_(eulerSolver), heatSolver_(heatSolver) { }

  void apply(AnyRange<CellIdx> eulerGhostCells) const override {
    for(auto eulerGhostCellId : eulerGhostCells) {
      CellIdx eulerBndryCellId, heatBndryCellId;
      std::tie(eulerBndryCellId, heatBndryCellId)
          = local_bndry_ids(eulerGhostCellId,eulerSolver_,heatSolver_);

      // velocity
      for(SInd d = 0; d < EulerSolver::nd; ++d) {
        eulerSolver_.cells().lhs(eulerGhostCellId,EulerSolver::V::rho_u(d))
            = 2.0 * 0.0 - eulerSolver_.cells().lhs(eulerBndryCellId,EulerSolver::V::rho_u(d));
        eulerSolver_.cells().rhs(eulerGhostCellId,EulerSolver::V::rho_u(d))
            = 2.0 * 0.0 - eulerSolver_.cells().rhs(eulerBndryCellId,EulerSolver::V::rho_u(d));
      }

      // density
      eulerSolver_.cells().lhs(eulerGhostCellId,EulerSolver::V::rho())
          = eulerSolver_.cells().lhs(eulerBndryCellId,EulerSolver::V::rho());
      eulerSolver_.cells().rhs(eulerGhostCellId,EulerSolver::V::rho())
          = eulerSolver_.cells().rhs(eulerBndryCellId,EulerSolver::V::rho());

      // energy
      const Num surfaceEnergy = 0.5 * (eulerSolver_.cells().lhs(eulerBndryCellId,EulerSolver::V::E())
                                       + heatSolver_.cells().lhs(heatBndryCellId,HeatSolver::V::T()));
      eulerSolver_.cells().lhs(eulerGhostCellId,EulerSolver::V::E())
          = 2 * surfaceEnergy  - eulerSolver_.cells().lhs(eulerBndryCellId,EulerSolver::V::E());
      eulerSolver_.cells().rhs(eulerGhostCellId,EulerSolver::V::E())
          = 2 * surfaceEnergy - eulerSolver_.cells().rhs(eulerBndryCellId,EulerSolver::V::E());
    }
  }

 private:
  EulerSolver& eulerSolver_;
  const HeatSolver& heatSolver_;
};


template<class HeatSolver>
struct HeatHeatWall : boundary::Condition {
  HeatHeatWall(HeatSolver& heatSolver0, HeatSolver& heatSolver1)
      : heatSolver0_(heatSolver0), heatSolver1_(heatSolver1) { }

  void apply(AnyRange<CellIdx> heat0GhostCells) const override {
    for(auto heat0GhostCellId : heat0GhostCells) {
      Ind heat0BndryCellId, localHeat1BndryCellId;
      std::tie(heat0BndryCellId, localHeat1BndryCellId)
          = local_bndry_ids(heat0GhostCellId,heatSolver0_,heatSolver1_);

      Num surfaceTemperature = 0.5 * (heatSolver0_.cells().lhs(heat0BndryCellId,0)
                                 + heatSolver1_.cells().lhs(localHeat1BndryCellId,0));
      heatSolver0_.cells().lhs(heat0GhostCellId,0)
          = 2.0 * surfaceTemperature - heatSolver0_.cells().lhs(heat0BndryCellId,0);
      heatSolver0_.cells().rhs(heat0GhostCellId,0)
          = 2.0 * surfaceTemperature - heatSolver0_.cells().rhs (heat0BndryCellId,0);
    }
  }

 private:
  HeatSolver& heatSolver0_;
  HeatSolver& heatSolver1_;
};

template<class EulerSolver>
struct EulerEulerWall : boundary::Condition {
  EulerEulerWall(EulerSolver& eulerSolver0, const EulerSolver& eulerSolver1)
      : eulerSolver0_(eulerSolver0), eulerSolver1_(eulerSolver1) { }

  void apply(AnyRange<CellIdx> euler0GhostCells) const override {
    for(auto euler0GhostCellId : euler0GhostCells) {

      CellIdx euler0BndryCellId, localEuler1BndryCellId;
      std::tie(euler0BndryCellId, localEuler1BndryCellId)
          = local_bndry_ids(euler0GhostCellId,eulerSolver0_,eulerSolver1_);

      for(SInd v = 0; v < EulerSolver::nvars; ++v) {
        eulerSolver0_.cells().lhs(euler0GhostCellId,v)
            = eulerSolver1_.cells().lhs(localEuler1BndryCellId,v);
         eulerSolver0_.cells().rhs(euler0GhostCellId,v)
             = eulerSolver1_.cells().rhs(localEuler1BndryCellId,v);
      }
    }
  }

 private:
  EulerSolver& eulerSolver0_;
  const EulerSolver& eulerSolver1_;
};

} // namespace bc

}}} // solver::fv::coupling namespace

////////////////////////////////////////////////////////////////////////////////
#undef ENABLE_DBG_
////////////////////////////////////////////////////////////////////////////////
#endif
