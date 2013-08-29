#ifndef HOM3_SOLVERS_FV_HEAT_BOUNDARY_CONDITIONS_HPP_
#define HOM3_SOLVERS_FV_HEAT_BOUNDARY_CONDITIONS_HPP_
////////////////////////////////////////////////////////////////////////////////
/// \file \brief Implements the Heat Equantion's boundary conditions.
////////////////////////////////////////////////////////////////////////////////
#include "solver/fv/boundary_condition.hpp"
////////////////////////////////////////////////////////////////////////////////
namespace hom3 { namespace solver { namespace fv { namespace heat {
////////////////////////////////////////////////////////////////////////////////

/// \brief Heat-equation boundary conditions
namespace bc {

/// \brief Dirichlet boundary condition for the temperature
template<class Solver> struct Dirichlet : fv::bc::Condition<Dirichlet<Solver>> {
  /// \brief Imposes a \p temperature distribution
  template<class F> Dirichlet(Solver& solver, F&& temperature) noexcept
    : T_srfc(temperature), s(solver) {}

  /// \brief Imposes a constant surface \p temperature
  Dirichlet(Solver& solver, const Num temperature) noexcept
    : T_srfc([=](CellIdx) { return temperature; }), s(solver) {}

  /// \brief Applies the boundary condition to a range of ghost cells
  template<class _>
  void operator()(_, Range<CellIdx> ghost_cells) const noexcept {
    for (auto ghostIdx : ghost_cells) {
      const CellIdx bndryIdx = s.boundary_info(ghostIdx).bndryIdx;
      s.T(_(), ghostIdx)
          = this->dirichlet(s.T(_(), bndryIdx), T_srfc(ghostIdx));
    }
  }

  /// Surface temperature distribution
  const std::function<Num(CellIdx)> T_srfc;
  Solver& s;
};

/// \brief Neumann boundary condition for the temperature
template<class Solver> struct Neumann : fv::bc::Condition<Neumann<Solver>> {
  /// \brief Imposes a \p heatFlux distribution
  template<class F> Neumann(Solver& solver, F&& heatFlux) noexcept
    : g_srfc([=](CellIdx c) { return heatFlux(c, 0); }), s(solver) {}

  /// \brief Imposes a constant surface \p heatFlux
  Neumann(Solver& solver, const Num heatFlux) noexcept
    : g_srfc([=](CellIdx) { return heatFlux; }), s(solver) {}

  /// \brief Applies the boundary condition to a range of ghost cells
  template<class _>
  void operator()(_, Range<CellIdx> ghost_cells) const noexcept {
    for (auto ghostIdx : ghost_cells) {
      const CellIdx bndryIdx = s.boundary_info(ghostIdx).bndryIdx;
      s.T(_(), ghostIdx)
          = this->neumann(s.T(_(), bndryIdx), g_srfc(bndryIdx),
                          s.cell_dx(bndryIdx, ghostIdx));
    }
  }

  /// Surface flux distribution
  const std::function<Num(CellIdx)> g_srfc;
  Solver& s;
};

} // namespace bc

////////////////////////////////////////////////////////////////////////////////
} // namespace heat
} // namespace fv
} // namespace solver
} // namespace hom3
////////////////////////////////////////////////////////////////////////////////
#endif
