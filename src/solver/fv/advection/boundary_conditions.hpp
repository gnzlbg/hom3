#ifndef HOM3_SOLVERS_FV_ADVECTION_BOUNDARY_CONDITIONS_HPP_
#define HOM3_SOLVERS_FV_ADVECTION_BOUNDARY_CONDITIONS_HPP_
////////////////////////////////////////////////////////////////////////////////
/// \file \brief Implements the Advection equantion's boundary conditions.
////////////////////////////////////////////////////////////////////////////////
#include "solver/fv/boundary_condition.hpp"
////////////////////////////////////////////////////////////////////////////////
namespace hom3 { namespace solver { namespace fv { namespace advection {
////////////////////////////////////////////////////////////////////////////////

/// \brief Advection equation boundary conditions
namespace bc {

/// \brief Dirichlet boundary condition
template<class Solver> struct Dirichlet : fv::bc::Condition<Dirichlet<Solver>> {
  /// \brief Imposes a surface \p distribution
  template<class F> Dirichlet(Solver& solver, F&& distribution) noexcept
    : q_srfc(distribution), s(solver) {}

  /// \brief Imposes a constant surface \p value
  Dirichlet(Solver& solver, const Num value) noexcept
      : Dirichlet(solver, [=](CellIdx) { return value; }) {}

  /// \brief Applies the boundary condition to a range of ghost cells
  template<class _>
  void operator()(_, Range<CellIdx> ghost_cells) const noexcept {
    for (auto ghostIdx : ghost_cells) {
      const CellIdx bndryIdx = s.boundary_info(ghostIdx).bndryIdx;
      s.q(_(), ghostIdx)
          = this->dirichlet(s.q(_(), bndryIdx), q_srfc(ghostIdx));
    }
  }

  /// Surface distribution
  const std::function<Num(CellIdx)> q_srfc;
  Solver& s;
};

/// \brief Neumann boundary condition
template<class Solver> struct Neumann : fv::bc::Condition<Neumann<Solver>> {
  /// \brief Imposes a \p flux distribution
  template<class F> Neumann(Solver& solver, F&& flux) noexcept
    : g_srfc(flux), s(solver) {}

  /// \brief Imposes a constant surface \p flux
  Neumann(Solver& solver, const Num flux) noexcept
    : Neumann(solver, [=](CellIdx) { return flux; }) {}

  /// \brief Applies the boundary condition to a range of ghost cells
  template<class _>
  void operator()(_, Range<CellIdx> ghost_cells) const noexcept {
    for (auto ghostIdx : ghost_cells) {
      const CellIdx bndryIdx = s.boundary_info(ghostIdx).bndryIdx;
      s.q(_(), ghostIdx)
          = this->neumann(s.q(_(), bndryIdx), g_srfc(ghostIdx),
                          s.cell_dx(bndryIdx, ghostIdx));
    }
  }

  /// Surface flux distribution
  const std::function<Num(CellIdx)> g_srfc;
  Solver& s;
};

/// \brief Dirichlet boundary condition
template<class Solver>
struct Characteristic : fv::bc::Condition<Characteristic<Solver>> {
  Characteristic(Solver& solver) noexcept : s(solver) {}

  /// \brief Applies the boundary condition to a range of ghost cells
  template<class _>
  void operator()(_, Range<CellIdx> ghost_cells) const noexcept {
    using namespace container::hierarchical;
    for (auto ghostIdx : ghost_cells) {
      const auto bInfo = s.boundary_info(ghostIdx);
      const CellIdx bndryIdx = bInfo.bndryIdx;

      const auto pos = bInfo.ghostPos;
      const auto dir = neighbor_direction(pos);
      const auto v = s.velocity(bndryIdx)(dir);
      if ((v > 0. && is_neighbor_at(pos, neg_dir))
          || (v < 0. && is_neighbor_at(pos, pos_dir))) {
      s.q(_(), ghostIdx)
          = this->dirichlet(s.q(_(), bndryIdx), s.q(_(), bndryIdx));
      } else {
        s.q(_(), ghostIdx)
            = this->neumann(s.q(_(), bndryIdx), 0.0,
                            s.cell_dx(bndryIdx, ghostIdx));
      }
    }
  }

  Solver& s;
};

}  // namespace bc

////////////////////////////////////////////////////////////////////////////////
}  // namespace advection
}  // namespace fv
}  // namespace solver
}  // namespace hom3
////////////////////////////////////////////////////////////////////////////////
#endif
