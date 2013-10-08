#ifndef HOM3_SOLVERS_FV_CNS_BOUNDARY_CONDITIONS_HPP_
#define HOM3_SOLVERS_FV_CNS_BOUNDARY_CONDITIONS_HPP_
////////////////////////////////////////////////////////////////////////////////
/// \file \brief Implements the Navier-Stokes boundary conditions.
////////////////////////////////////////////////////////////////////////////////
#include "solver/fv/boundary_condition.hpp"
#include "initial_conditions.hpp"
#include "solver/fv/euler/boundary_conditions.hpp"
////////////////////////////////////////////////////////////////////////////////
namespace hom3 { namespace solver { namespace fv { namespace cns {
////////////////////////////////////////////////////////////////////////////////

/// \brief Navier-Stokes equations boundary conditions
namespace bc {

/// \brief Neumann boundary condition (all gradients)
template<class Solver> struct Neumann : euler::bc::Neumann<Solver> {
  using euler::bc::Neumann<Solver>::template Neumann;
};

/// \brief Boundary conditions for external flows
namespace external {

/// \brief Inflow boundary condition
///
/// Prescribes $u_\Gamma = u_D$, $\rho_\Gamma = \rho_D$,
/// $\nabla_n p \vert_\Gamma = 0$.
template<class Solver> struct Inflow : fv::bc::Condition<Inflow<Solver>> {
  static const SInd nd = Solver::nd;
  static const SInd nvars = Solver::nvars;
  using V = Indices<nd>;

  /// Distribution returns $u_D$ and $\rho_D$
  template<class F> Inflow(Solver& solver, F&& distribution) noexcept
    : condition(distribution), s(solver) {}

  template<class _>
  void operator()(_, Range<CellIdx> ghost_cells) const noexcept {
    for (auto ghostIdx : ghost_cells) {
      const CellIdx bndryIdx = s.boundary_info(ghostIdx).bndryIdx;
      const auto bndryPVars = s.pv(_(), bndryIdx);  /// PV of boundary cell

      /// Surface primitive variables
      const Num rhoD = condition(bndryIdx, V::rho());
      NumA<nd> uD;
      for (auto d : s.grid().dimensions()) {
        uD(d) = condition(bndryIdx, V::u(d));
      }

      NumA<nvars> gCellVars;  /// PV of ghost cell:
      for (auto d : s.grid().dimensions()) {  /// velocity_surface = uD
        gCellVars(V::u(d)) = this->dirichlet(bndryPVars(V::u(d)), uD(d));
      }
      /// density_surface = rhoD
      gCellVars(V::rho()) = this->dirichlet(bndryPVars(V::rho()), rhoD);
      /// normal pressure gradient = 0
      gCellVars(V::p()) = this->neumann(bndryPVars(V::p()));

      s.Q(_(), ghostIdx) = s.cv(gCellVars);
    }
  }

  std::function<Num(CellIdx, SInd)> condition;
  Solver& s;
};


/// \brief Outflow boundary condition
///
/// Prescribes $p_\gamma = p_D$, $\nabla_n \mathbf{u} \vert_\Gamma = 0$,
/// $\nabla_n \rho \vert_\Gamma = 0$.
template<class Solver> struct Outflow : fv::bc::Condition<Outflow<Solver>>  {
  static const SInd nd = Solver::nd;
  static const SInd nvars = Solver::nvars;
  using V = Indices<nd>;

  /// \brief Distribution returns $p_D$
  template<class F> Outflow(Solver& solver, F&& distribution) noexcept
    : condition(distribution), s(solver) {}

  template<class _>
  void operator()(_, Range<CellIdx> ghost_cells) const noexcept {
    for (auto ghostIdx : ghost_cells) {
      const CellIdx bndryIdx = s.boundary_info(ghostIdx).bndryIdx;
      const auto bndryPVars = s.pv(_(), bndryIdx);  /// PV of boundary cell

      const Num pD = condition(bndryIdx, V::p());

      NumA<nvars> gCellVars;  /// PV of ghost cell:

      for (auto d : s.grid().dimensions()) {  // normal velocity gradient = 0
        gCellVars(V::u(d)) = this->neumann(bndryPVars(V::u(d)));
      }
      /// normal density gradient = 0
      gCellVars(V::rho()) = this->neumann(bndryPVars(V::rho()));
      /// pressure_surface = pD
      gCellVars(V::p()) = this->dirichlet(bndryPVars(V::p()), pD);

      s.Q(_(), ghostIdx) = s.cv(gCellVars);
    }
  }

  std::function<Num(CellIdx, SInd)> condition;
  Solver& s;
};

}  // namespace external

/// \brief Wall boundary conditions
namespace wall {

/// \brief Adiabatic no-slip wall
///
/// Prescribes $\nabla_n p \vert_\Gamma = \nabla_n \rho \vert_\Gamma = 0$,
/// $\mathbf{u}_\Gamma = 0$.
template<class Solver>
struct AdiabaticNoSlip : fv::bc::Condition<AdiabaticNoSlip<Solver>>  {
  static const SInd nd = Solver::nd;
  static const SInd nvars = Solver::nvars;
  using V = Indices<nd>;

  explicit AdiabaticNoSlip(Solver& solver) : s(solver) {}

  template<class _>
  void operator()(_, Range<CellIdx> ghost_cells) const noexcept {
    for (auto ghostIdx : ghost_cells) {
      const CellIdx bndryIdx = s.boundary_info(ghostIdx).bndryIdx;
      const auto bndryPVars = s.pv(_(), bndryIdx);  /// PV of boundary cell

      NumA<nvars> gCellVars;  /// PV of ghost cell:
      for (auto d : s.grid().dimensions()) {  /// velocity_srfc = 0
        gCellVars(V::u(d)) = this->dirichlet(bndryPVars(V::u(d)));
      }
      /// density gradient normal to surface = 0
      gCellVars(V::rho()) = this->neumann(bndryPVars(V::rho()));
      /// pressure gradient normal to surface = 0
      gCellVars(V::p()) = this->neumann(bndryPVars(V::p()));

      s.Q(_(), ghostIdx) = s.cv(gCellVars);
    }
  }

  template<class _> Num slope
  (const CellIdx bndryIdx, const SInd v, const SInd dir) const noexcept {
    const auto tmp = s.template slope<_>(bndryIdx, v, dir);
    return v < nd ? - tmp : tmp;  // zero tangential gradient for velocity
  }

  Solver& s;
};

/// \brief Adiabatic slip wall
///
/// Prescribes $\nabla_n p \vert_\Gamma = \nabla_n \rho \vert_\Gamma = \nabla_t
/// \mathbf{u} \vert_\Gamma = 0$, $\mathbf{u} \cdot \mathbf{n} = 0$.
template<class Solver>
struct AdiabaticSlip  : fv::bc::Condition<AdiabaticSlip<Solver>> {
  static const SInd nd = Solver::nd;
  static const SInd nvars = Solver::nvars;
  using V = Indices<nd>;

  explicit AdiabaticSlip(Solver& solver) : s(solver) {}

  template<class _>
  void operator()(_, Range<CellIdx> ghost_cells) const noexcept {
    for (auto ghostIdx : ghost_cells) {
      const auto bndryInfo = s.boundary_info(ghostIdx);
      const CellIdx bndryIdx = bndryInfo.bndryIdx;
      auto pos_wrt_ghostCell = bndryInfo.bndryPos;

      auto normalDirection
        = static_cast<SInd>(std::floor(pos_wrt_ghostCell / 2));
      const auto bndryPVars = s.pv(_(), bndryIdx);  /// PV of boundary cell

      NumA<nvars> gCellVars;  /// PV of ghost cell:

      // Neumann BC: velocity
      for (auto d : s.grid().dimensions()) {  /// normal velocity gradient = 0
        gCellVars(d) = this->neumann(bndryPVars(d));
      }
      /// Except for the normal component: normal velocity_srfc = 0
      gCellVars(normalDirection) = this->dirichlet(bndryPVars(normalDirection));
      // normal density gradient = 0
      gCellVars(V::rho()) = this->neumann(bndryPVars(V::rho()));
      // normal pressure gradient = 0
      gCellVars(V::p()) = this->neumann(bndryPVars(V::p()));

      s.Q(_(), ghostIdx) = s.cv(gCellVars);
    }
  }

  Solver& s;
};

/// \brief Isothermal no-slip wall
///
/// Prescribes $\nabla_n p \vert_\Gamma = \nabla_n \rho \vert_\Gamma = 0$,
/// $\mathbf{u}_\Gamma = 0$.
template<class Solver>
struct IsothermalNoSlip : fv::bc::Condition<IsothermalNoSlip<Solver>> {
  static const SInd nd = Solver::nd;
  static const SInd nvars = Solver::nvars;
  using V = Indices<nd>;

  /// \brief Temperature $T_D$
  template<class F> IsothermalNoSlip(Solver& solver, F&& temperature) noexcept
    : T_srfc(temperature), s(solver) {}

  template<class _>
  void operator()(_, Range<CellIdx> ghost_cells) const noexcept {
    for (auto ghostIdx : ghost_cells) {
      const CellIdx bndryIdx = s.boundary_info(ghostIdx).bndryIdx;
      const auto bndryPVars = s.pv(_(), bndryIdx);  /// PV of boundary cell

      NumA<nvars> gCellVars;  /// PV of ghost cell:
      for (auto d : s.grid().dimensions()) {  /// velocity surface = 0
        gCellVars(d) = this->dirichlet(bndryPVars(d));
      }
      /// density surface = gamma * p_srfc / T_srfc
      auto p_srfc = bndryPVars(V::p());  // i.e. ~= p_bndryCell
      const Num rho_srfc = p_srfc * s.quantities.gamma() / T_srfc(ghostIdx);
      gCellVars(V::rho()) = this->dirichlet(bndryPVars(V::rho()), rho_srfc);
      /// normal pressure gradient = 0
      gCellVars(V::p()) = this->neumann(bndryPVars(V::p()));

      s.Q(_(), ghostIdx) = s.cv(gCellVars);
    }
  }

  template<class _> Num slope
  (const CellIdx bndryIdx, const SInd v, const SInd dir) const noexcept {
    const auto tmp = s.template slope<_>(bndryIdx, v, dir);
    return v < nd ? - tmp : tmp;  // zero tangential gradient for velocity
  }

  std::function<Num(CellIdx)> T_srfc;
  Solver& s;
};

}  // namespace wall

}  // namespace bc

////////////////////////////////////////////////////////////////////////////////
}  // namespace cns
}  // namespace fv
}  // namespace solver
}  // namespace hom3
////////////////////////////////////////////////////////////////////////////////
#endif
