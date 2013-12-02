#ifndef HOM3_GRID_BOUNDARY_HPP_
#define HOM3_GRID_BOUNDARY_HPP_
////////////////////////////////////////////////////////////////////////////////
/// \file \brief Grid boundary interface
////////////////////////////////////////////////////////////////////////////////
/// Includes:
#include <memory>
#include <functional>
#include "globals.hpp"
/// Options:
#define ENABLE_DBG_ 0
#include "misc/dbg.hpp"
////////////////////////////////////////////////////////////////////////////////
namespace hom3 { namespace grid {
////////////////////////////////////////////////////////////////////////////////

/// \brief Boundary related grid functionality
namespace boundary {

/// \brief Implements the grid boundary concept
///
/// Each grid boundary has
///   - an associated geometry,
/// and allows
///   - queriying the signed-distance field of the boundary.
///
/// Note: the geometry shared_ptr is capture by value in the lambda
template<SInd nd> struct Interface {
  /// \brief Constructs a boundary with a custom boundary condition
  template<class Geometry, class Solver>
  Interface(const String name,
            std::shared_ptr<Geometry> geometry,
            const Solver& solver) noexcept
    : signed_distance([=](const NumA<nd> x) {
        return (*geometry)(x);
      })
    , name_(name)
    , solverIdx_(solver.solver_idx())
  {}

  /// Signed-distance to the boundary
  const std::function<Num(const NumA<nd>&)> signed_distance;
  /// \brief Index of the solver owning the boundary condition
  inline SolverIdx solver_idx() const noexcept { return solverIdx_; }
  /// \brief Boundary condition name
  inline String name() const noexcept { return name_; }
 private:
  /// Boundary name
  const String name_;
  /// Index of the solver owning the boundary condition
  const SolverIdx solverIdx_;
};

}  // namespace boundary

////////////////////////////////////////////////////////////////////////////////
}  // namespace grid
}  // namespace hom3
////////////////////////////////////////////////////////////////////////////////
#undef ENABLE_DBG_
////////////////////////////////////////////////////////////////////////////////
#endif
