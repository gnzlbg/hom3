#ifndef HOM3_GEOMETRY_BOUNDARY_HPP_
#define HOM3_GEOMETRY_BOUNDARY_HPP_
////////////////////////////////////////////////////////////////////////////////
#include <string>
#include <vector>
#include "../globals.hpp"
#include "../solver/interface.hpp"
////////////////////////////////////////////////////////////////////////////////
/// \file \brief This file implements the boundary condition interface
////////////////////////////////////////////////////////////////////////////////

namespace hom3 {

/// \brief Boundary conditions
namespace boundary {

/// \brief Boundary condition interface
///
/// Single solver BC: apply bc to solver's cells
///
/// Solver/Solver BC: apply bc to active_solver's cells using information from
/// passive_solver
struct Condition {
  void operator()(AnyRange<CellIdx> cells) const { apply(cells); }
 private:
  virtual void apply(AnyRange<CellIdx> /*cells*/) const {
    ASSERT(false, "calling an empty bc!");
  };
};

using ConditionHandler = std::shared_ptr<Condition>;

template<class T, class... Args>
ConditionHandler make_condition(Args&&... args) {
  return std::make_shared<T>(std::forward<Args>(args)...);
}

/// \brief Implements the boundary concept
///
/// Each boundary has:
/// - an associated geometry, and
/// - a boundary condition functor.
///
/// Each boundary allows:
/// - querying the signed-distance field of the boundary, and
/// - applying its boundary condition for a given cell of a given solver, by
/// either:
///   - using the solver interface (calling a standard boundary condition), or
///   - using a user-provided boundary condition functor.
///
/// Note: the geometry shared_ptr is capture by value in the lambda which is
/// then stored in the std::function, i.e. this gives value semantics to the
/// boundary Interface (even when all other references to the geometry have been
/// destroyed, the interface still owns the object :)
///
/// \todo rename to hom3::Boundary (it's not an interface, is a boundary with
/// value semantics)
/// \todo change shared_ptr to unique_ptr for the geometry
/// object geometry when C++14 (requires lambda capture-by-move)
template<SInd nd> struct Interface {
  using Conditions = std::vector<ConditionHandler>;
  /// \brief Constructs a boundary with a custom boundary condition
  template<class Geometry>
  Interface(const std::string bcName, std::shared_ptr<Geometry> geometry,
            const SolverIdx solverIdx,
            const std::vector<ConditionHandler> bc)
      : name_(bcName),
        signed_distance([=](const NumA<nd> x){ return geometry->operator()(x); }),
        conditions_(bc),
        solverIdx_(solverIdx)
  {}

  Interface() = default;
  Interface(const Interface&) = default;
  Interface& operator=(const Interface&) = default;

  std::string name_;
  /// \brief Signed-distance to the boundary
  std::function<Num(const NumA<nd>&)> signed_distance;
  /// \brief Applies boundary condition on boundary cell \p localBCellId of
  /// solver \p s
  inline std::string name() const { return name_; }
  inline SolverIdx solver_idx() const { return solverIdx_; }
  inline bool is_valid(const SolverIdx solverIdx) const {return solverIdx == solver_idx();}
  inline const Condition& condition(const SInd condition) const {return (*conditions_[condition]);}
  inline const Conditions& conditions() const { return conditions_; }
 private:
  Conditions conditions_;
  const SolverIdx solverIdx_;
};

}} // hom3::boundary namespace

////////////////////////////////////////////////////////////////////////////////
#endif
