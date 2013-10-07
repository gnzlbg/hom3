#ifndef HOM3_SOLVERS_FV_EULER_BOUNDARY_CONDITIONS_HPP_
#define HOM3_SOLVERS_FV_EULER_BOUNDARY_CONDITIONS_HPP_
////////////////////////////////////////////////////////////////////////////////
/// \file \brief Implements the Euler-equations' boundary conditions.
////////////////////////////////////////////////////////////////////////////////
#include "solver/fv/boundary_condition.hpp"
#include "initial_conditions.hpp"
////////////////////////////////////////////////////////////////////////////////
namespace hom3 { namespace solver { namespace fv { namespace euler {
////////////////////////////////////////////////////////////////////////////////

/// \brief Euler-equations boundary conditions
namespace bc {

/// \brief Neumann boundary condition for all variables:
/// $\nabla_n \mathbf{Q} \vert_\Gamma = \mathbf{g}
template<class Solver> struct Neumann : fv::bc::Condition<Neumann<Solver>> {
  /// \brief Imposes a \p flux distribution
  template<class F> Neumann(Solver& solver, F&& flux) noexcept
    : g_srfc(flux), s(solver) {}

  /// \brief Imposes a constant surface \p flux
  Neumann(Solver& solver, const NumA<Solver::nvars> flux
                          = NumA<Solver::nvars>::Zero())
    : g_srfc([=](CellIdx, const SInd v) { return flux(v); }), s(solver) {}

  /// \brief Applies the boundary condition to a range of ghost cells
  template<class _>
  void operator()(_, Range<CellIdx> ghost_cells) const noexcept {
    for (auto ghostIdx : ghost_cells) {
      const CellIdx bndryIdx = s.boundary_info(ghostIdx).bndryIdx;
      const auto dx = s.cell_dx(bndryIdx, ghostIdx);
      for (SInd v = 0; v < s.nvars; ++v) {
        s.Q(_(), ghostIdx, v)
          = this->neumann(s.Q(_(), bndryIdx, v), g_srfc(ghostIdx, v), dx);
      }
    }
  }

  /// Surface flux distribution
  std::function<Num(CellIdx, SInd)> g_srfc;
  Solver& s;
};

/// \brief Dirichlet boundary condition for all variables:
/// $\mathbf{Q} \vert_\Gamma = \mathbf{Q} \vert_D
template<class Solver> struct Dirichlet : fv::bc::Condition<Dirichlet<Solver>> {
  /// \brief Imposes a \p flux distribution
  template<class F> Dirichlet(Solver& solver, F&& variables) noexcept
    : Q_srfc(variables), s(solver) {}

  /// \brief Applies the boundary condition to a range of ghost cells
  template<class _>
  void operator()(_, Range<CellIdx> ghost_cells) const noexcept {
    for (auto ghostIdx : ghost_cells) {
      const CellIdx bndryIdx = s.boundary_info(ghostIdx).bndryIdx;
      for (SInd v = 0; v < s.nvars; ++v) {
        s.Q(_(), ghostIdx, v)
          = this->dirichlet(s.Q(_(), bndryIdx, v), Q_srfc(ghostIdx, v));
      }
    }
  }

  /// Surface variables
  std::function<Num(CellIdx, SInd)> Q_srfc;
  Solver& s;
};

/// \brief Analytical boundary condition for the isentropic vortex case
template<class Solver>
struct IsentropicVortex : fv::bc::Condition<IsentropicVortex<Solver>> {
  static const auto nd = Solver::nd;
  static const auto nvars = Solver::nvars;
  using V = Indices<Solver::nd>;
  explicit IsentropicVortex(Solver& solver) noexcept : s(solver) {}

  /// \brief Applies the boundary condition to a range of ghost cells
  template<class _>
  void operator()(_, Range<CellIdx> ghost_cells) const noexcept {
    for (auto ghostIdx : ghost_cells) {
      const CellIdx bndryIdx = s.boundary_info(ghostIdx).bndryIdx;

      const NumA<nd>    x_bc = s.cells().x_center.row(bndryIdx);
      const NumA<nd>    x_gc = s.cells().x_center.row(ghostIdx);
      const NumA<nvars> gcVars = ic::isentropic_vortex<nd>(x_gc, s.time());
      const NumA<nvars> bcVars = ic::isentropic_vortex<nd>(x_bc, s.time());
      /// \todo vars = vars
      for (SInd v = 0; v < nvars; ++v) {
        s.Q(_(), ghostIdx, v) = gcVars(v);
        s.Q(_(), bndryIdx, v) = bcVars(v);
      }
    }
  }

  Solver& s;
};

}  // namespace bc

////////////////////////////////////////////////////////////////////////////////
}  // namespace euler
}  // namespace fv
}  // namespace solver
}  // namespace hom3
////////////////////////////////////////////////////////////////////////////////
#endif
