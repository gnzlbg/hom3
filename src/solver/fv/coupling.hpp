#ifndef HOM3_SOLVERS_FV_COUPLING_HPP_
#define HOM3_SOLVERS_FV_COUPLING_HPP_
////////////////////////////////////////////////////////////////////////////////
/// \file \brief This file implements the solver coupling functionality.
////////////////////////////////////////////////////////////////////////////////
/// Includes:
#include <algorithm>
#include <type_traits>
#include "solver/fv/heat/tags.hpp"
#include "solver/fv/cns/tags.hpp"
////////////////////////////////////////////////////////////////////////////////
/// Options:
#define ENABLE_DBG_ 0
#include "misc/dbg.hpp"
////////////////////////////////////////////////////////////////////////////////

namespace hom3 { namespace solver { namespace fv {

namespace coupling {

/// \brief Computes the local boundary cell ids of solver0 and solver1
///
/// \returns (solver0LocalBdryCellId, solver1LocalBndryCellId)
template<class Solver0, class Solver1>
std::tuple<CellIdx, CellIdx> local_bndry_ids
(const CellIdx ghostIdx0, const Solver0& s0, const Solver1& s1) {
  ASSERT(is_valid(ghostIdx0), "invalid ghost cell idx!");
  const auto bndryInfo0 = s0.boundary_info(ghostIdx0);
  const auto bndryIdx0 = bndryInfo0.bndryIdx;
  const auto ghostPos0 = bndryInfo0.ghostPos;

  const auto bndryNodeIdx0 = s0.node_idx(bndryIdx0);
  const auto bndryNodeIdx1 = s0.grid().find_samelvl_neighbor
                             (bndryNodeIdx0, ghostPos0);
  const auto bndryIdx1 = s0.grid().cell_idx(bndryNodeIdx1, s1.solver_idx());
  ASSERT(is_valid(bndryIdx0), "invalid boundary cell idx!");
  ASSERT(is_valid(bndryIdx1), "invalid boundary cell idx!");
  return std::make_tuple(bndryIdx0, bndryIdx1);
}

/// \brief Couples a Navier-Stokes solver for ideal compressible flow with a
/// Heat equation solver.
///
/// $\frac{ t_\mathrm{solid} }{ t_\mathrm{fluid} } = \frac{
/// \overline{L}_\mathrm{ref,solid}^2 }{ \overline{L}_\mathrm{ref,fluid}^2 }
/// \mathrm{Tr} \mathrm{Re}_0 \mathrm{Pr}_0 = \mathrm{scaling} $
template<class CNSSolver, class HeatSolver,
  EnableIf<std::is_same<typename CNSSolver::physics_type,
                        cns::type_tag>>  = traits::dummy,
  EnableIf<std::is_same<typename HeatSolver::physics_type,
                        heat::type_tag>> = traits::dummy>
inline void apply
(const Num Tr0, CNSSolver& cnsSolver, HeatSolver& heatSolver) noexcept {
  const Num scaling = Tr0 * cnsSolver.quantities.Pr0()
                      * cnsSolver.quantities.Re0();

  const Num min_dt_cns = cnsSolver.min_dt();
  const Num min_dt_heat = heatSolver.min_dt() / scaling;
  const Num dt = std::min(min_dt_cns, min_dt_heat);

  cnsSolver.dt(dt);
  heatSolver.dt(dt * scaling);
}

/// \brief Couples two equal solvers by taking the minimum time-step from both.
template<class Solver> inline void apply(Solver& s1, Solver& s2) noexcept {
  const Num min_dt_s1 = s1.min_dt();
  const Num min_dt_s2 = s2.min_dt();
  const Num dt = std::min(min_dt_s1, min_dt_s2);

  s1.dt(dt);
  s2.dt(dt);
}

}  // namespace coupling

////////////////////////////////////////////////////////////////////////////////
}  // namespace fv
}  // namespace solver
}  // namespace hom3
////////////////////////////////////////////////////////////////////////////////
#undef ENABLE_DBG_
////////////////////////////////////////////////////////////////////////////////
#endif
