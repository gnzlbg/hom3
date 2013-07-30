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
std::tuple<Ind,Ind> local_bndry_ids(const Ind solver0GhostCellId,
                                    Solver0& solver0, Solver1& solver1) {

  Ind solver0BndryCellId; SInd solver0BndryCellPos;
  std::tie(solver0BndryCellId,solver0BndryCellPos)
      = solver0.boundary_cell_id(solver0GhostCellId);

  SInd solver0GhostCellPos
      = solver0.grid().nodes().opposite_nghbr_position(solver0BndryCellPos);

  Ind solver0GlobalBndryCellId = solver0.global_id(solver0BndryCellId);
  Ind solver1GlobalBndryCellId
      = solver0.grid().nodes().find_samelvl_nghbr(solver0GlobalBndryCellId,
                                                        solver0GhostCellPos);
  Ind solver1LocalBndryCellId
      = solver0.grid().nodes().local_id(solver1GlobalBndryCellId,
                                             solver1.solver_id());
  return std::make_tuple(solver0BndryCellId, solver1LocalBndryCellId);
}

namespace bc {

template<class HeatSolver, class EulerSolver>
struct HeatEulerWall : boundary::Condition {
  HeatEulerWall(HeatSolver& heatSolver, EulerSolver& eulerSolver)
      : heatSolver_(heatSolver), eulerSolver_(eulerSolver) { }

  void apply(AnyRange<Ind> heatGhostCells) const override {
    for(auto heatGhostCellId : heatGhostCells) {
      Ind heatBndryCellId, localEulerBndryCellId;
      std::tie(heatBndryCellId, localEulerBndryCellId)
          = local_bndry_ids(heatGhostCellId,heatSolver_,eulerSolver_);

      Num surfaceEnergy = 0.5 * (heatSolver_.cells().vars(heatBndryCellId,0)
                                 + eulerSolver_.cells().vars(localEulerBndryCellId,
                                                             EulerSolver::V::E()));
      heatSolver_.cells().vars(heatGhostCellId,0)
          = 2.0 * surfaceEnergy - heatSolver_.cells().vars(heatBndryCellId,0);
      heatSolver_.cells().rhs(heatGhostCellId,0)
          = 2.0 * surfaceEnergy - heatSolver_.cells().rhs (heatBndryCellId,0);
    }
  }

 private:
  HeatSolver& heatSolver_;
  EulerSolver& eulerSolver_;
};


template<class EulerSolver, class HeatSolver>
struct EulerHeatWall : boundary::Condition {

  EulerHeatWall(EulerSolver& eulerSolver, HeatSolver& heatSolver)
      : eulerSolver_(eulerSolver), heatSolver_(heatSolver) { }

  void apply(AnyRange<Ind> eulerGhostCells) const override {
    for(auto eulerGhostCellId : eulerGhostCells) {
      Ind eulerBndryCellId, localHeatBndryCellId;
      std::tie(eulerBndryCellId, localHeatBndryCellId)
          = local_bndry_ids(eulerGhostCellId,eulerSolver_,heatSolver_);

      // velocity
      for(SInd d = 0; d < EulerSolver::nd; ++d) {
        eulerSolver_.cells().vars(eulerGhostCellId,EulerSolver::V::rho_u(d))
            = 2.0 * 0.0 - eulerSolver_.cells().vars(eulerBndryCellId,EulerSolver::V::rho_u(d));
        eulerSolver_.cells().rhs(eulerGhostCellId,EulerSolver::V::rho_u(d))
            = 2.0 * 0.0 - eulerSolver_.cells().rhs (eulerBndryCellId,EulerSolver::V::rho_u(d));
      }

      // density
      eulerSolver_.cells().vars(eulerGhostCellId,EulerSolver::V::rho())
          = eulerSolver_.cells().vars(eulerBndryCellId,EulerSolver::V::rho());
      eulerSolver_.cells().rhs(eulerGhostCellId,EulerSolver::V::rho())
          = eulerSolver_.cells().rhs (eulerBndryCellId,EulerSolver::V::rho());

      // energy
      Num surfaceEnergy = 0.5 * (eulerSolver_.cells().vars(eulerBndryCellId,EulerSolver::V::E())
                              + heatSolver_.cells().vars(localHeatBndryCellId));
      eulerSolver_.cells().vars(eulerGhostCellId,EulerSolver::V::E())
          = 2 * surfaceEnergy  - eulerSolver_.cells().vars(eulerBndryCellId,EulerSolver::V::E());
      eulerSolver_.cells().rhs(eulerGhostCellId,EulerSolver::V::E())
          = 2 * surfaceEnergy - eulerSolver_.cells().rhs (eulerBndryCellId,EulerSolver::V::E());
    }
  }

 private:
  EulerSolver& eulerSolver_;
  HeatSolver& heatSolver_;
};


template<class HeatSolver>
struct HeatHeatWall : boundary::Condition {
  HeatHeatWall(HeatSolver& heatSolver0, HeatSolver& heatSolver1)
      : heatSolver0_(heatSolver0), heatSolver1_(heatSolver1) { }

  void apply(AnyRange<Ind> heat0GhostCells) const override {
    for(auto heat0GhostCellId : heat0GhostCells) {
      Ind heat0BndryCellId, localHeat1BndryCellId;
      std::tie(heat0BndryCellId, localHeat1BndryCellId)
          = local_bndry_ids(heat0GhostCellId,heatSolver0_,heatSolver1_);

      Num surfaceTemperature = 0.5 * (heatSolver0_.cells().vars(heat0BndryCellId,0)
                                 + heatSolver1_.cells().vars(localHeat1BndryCellId,0));
      heatSolver0_.cells().vars(heat0GhostCellId,0)
          = 2.0 * surfaceTemperature - heatSolver0_.cells().vars(heat0BndryCellId,0);
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
  EulerEulerWall(EulerSolver& eulerSolver0, EulerSolver& eulerSolver1)
      : eulerSolver0_(eulerSolver0), eulerSolver1_(eulerSolver1) { }

  void apply(AnyRange<Ind> euler0GhostCells) const override {
    for(auto euler0GhostCellId : euler0GhostCells) {

      Ind euler0BndryCellId, localEuler1BndryCellId;
      std::tie(euler0BndryCellId, localEuler1BndryCellId)
          = local_bndry_ids(euler0GhostCellId,eulerSolver0_,eulerSolver1_);

      for(SInd v = 0; v < EulerSolver::nvars; ++v) {
        eulerSolver0_.cells().vars(euler0GhostCellId,v)
            = eulerSolver1_.cells().vars(localEuler1BndryCellId,v);
         eulerSolver0_.cells().rhs(euler0GhostCellId,v)
             = eulerSolver1_.cells().rhs(localEuler1BndryCellId,v);
      }
    }
  }

 private:
  EulerSolver& eulerSolver0_;
  EulerSolver& eulerSolver1_;
};

} // namespace bc

}}} // solver::fv::coupling namespace

////////////////////////////////////////////////////////////////////////////////
#undef ENABLE_DBG_
////////////////////////////////////////////////////////////////////////////////
#endif
