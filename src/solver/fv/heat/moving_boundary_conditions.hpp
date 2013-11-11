#ifndef HOM3_SOLVERS_FV_HEAT_MOVING_BOUNDARY_CONDITIONS_HPP_
#define HOM3_SOLVERS_FV_HEAT_MOVING_BOUNDARY_CONDITIONS_HPP_
////////////////////////////////////////////////////////////////////////////////
/// \file \brief Implements the Heat Equantion's moving boundary conditions.
////////////////////////////////////////////////////////////////////////////////
#include "solver/fv/boundary_condition.hpp"
////////////////////////////////////////////////////////////////////////////////
namespace hom3 { namespace solver { namespace fv { namespace heat {
////////////////////////////////////////////////////////////////////////////////
namespace bc {
////////////////////////////////////////////////////////////////////////////////

/// \brief Moving boundary conditions
namespace mb {

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
      s.T(_(), ghostIdx)
          = T_srfc(ghostIdx);
      // interpolate_values_to_moving_boundary_ghost_cell
      // (_(), ghostIdx, T_srfc(ghostIdx), 0);
    }
  }

  bool is_moving_boundary() const noexcept { return true; }

  /// Surface temperature distribution
  const std::function<Num(CellIdx)> T_srfc;
  Solver& s;
};


}  // namespace mb

////////////////////////////////////////////////////////////////////////////////
}  // namespace bc
}  // namespace heat
}  // namespace fv
}  // namespace solver
}  // namespace hom3
////////////////////////////////////////////////////////////////////////////////
#endif
