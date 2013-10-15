#ifndef HOM3_SOLVERS_FV_BOUNDARY_CONDITION_HPP_
#define HOM3_SOLVERS_FV_BOUNDARY_CONDITION_HPP_
////////////////////////////////////////////////////////////////////////////////
/// Includes:
#include "globals.hpp"
#include "grid/boundary.hpp"
#include "tags.hpp"
/// Options:
#define ENABLE_DBG_ 0
#include "misc/dbg.hpp"
////////////////////////////////////////////////////////////////////////////////
namespace hom3 { namespace solver { namespace fv {
////////////////////////////////////////////////////////////////////////////////

namespace bc {

/// \brief Implements the Finite Volume boundary condition concept
///
/// Extends the grid::boundary concept with
/// - an apply(GhostCellRange) function, that imposes the boundary condition
///
template<SInd nd> struct Interface : grid::boundary::Interface<nd> {
  using GhostCellRange = Range<CellIdx>;

  /// \brief Construct a boundary condition
  template<class Solver, class Geometry, class Condition>
  Interface(const String name, std::shared_ptr<Geometry> geometry,
            const Solver& solver, Condition&& condition) noexcept
      : grid::boundary::Interface<nd>(name, geometry, solver)
      , apply_lhs(condition), apply_rhs(condition)
      , slope_lhs([=](const CellIdx cIdx, const SInd v, const SInd dir) {
          return condition.template slope<lhs_tag>(cIdx, v, dir); })
      , slope_rhs([=](const CellIdx cIdx, const SInd v, const SInd dir) {
          return condition.template slope<rhs_tag>(cIdx, v, dir); })
  {}

  /// Applies the boundary condition to the lhs of a GhostCellRange
  void apply(lhs_tag, GhostCellRange ghost_cells) const noexcept
  { apply_lhs(lhs, ghost_cells); }

  /// Applies the boundary condition to the rhs of a GhostCellRange
  void apply(rhs_tag, GhostCellRange ghost_cells) const noexcept
  { apply_rhs(rhs, ghost_cells); }

  Num slope(lhs_tag, const CellIdx cIdx, const SInd v,
            const SInd dir) const noexcept
  { return slope_lhs(cIdx, v, dir); }

  Num slope(rhs_tag, const CellIdx cIdx, const SInd v,
            const SInd dir) const noexcept
  { return slope_rhs(cIdx, v, dir); }

 private:
  const std::function<void(lhs_tag, GhostCellRange)> apply_lhs;
  const std::function<void(rhs_tag, GhostCellRange)> apply_rhs;
  const std::function<Num(CellIdx, SInd, SInd)> slope_lhs;
  const std::function<Num(CellIdx, SInd, SInd)> slope_rhs;
};

/// \brief Contains general boundary condition functionality
template<class BC> struct Condition {
  /// \brief Imposes boundary cell's tangential slope on the boundary surface
  template<class _>
  Num slope(const CellIdx bndryIdx, const SInd v,
            const SInd dir) const noexcept {
    return static_cast<const BC*>(this)->s.template slope<_>(bndryIdx, v, dir);
  }

  /// \brief Ghost cell value for dirichlet boundary condition
  inline Num dirichlet(const Num boundaryCellValue,
                       const Num surfaceValue = 0) const noexcept
  { return 2.0 * surfaceValue - boundaryCellValue; }

  /// \brief Ghost cell value for dirichlet boundary condition
  inline Num neumann(const Num boundaryCellValue, const Num surfaceFlux = 0,
                     const Num cell_distance = 0) const noexcept
  { return boundaryCellValue - cell_distance * surfaceFlux; }
};

}  // namespace bc

////////////////////////////////////////////////////////////////////////////////
}  // namespace fv
}  // namespace solver
}  // namespace hom3
////////////////////////////////////////////////////////////////////////////////
#undef ENABLE_DBG_
////////////////////////////////////////////////////////////////////////////////
#endif
