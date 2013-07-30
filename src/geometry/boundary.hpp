#ifndef HOM3_GEOMETRY_BOUNDARY_HPP_
#define HOM3_GEOMETRY_BOUNDARY_HPP_
////////////////////////////////////////////////////////////////////////////////
#include "../globals.hpp"
#include "../solver/interface.hpp"
////////////////////////////////////////////////////////////////////////////////
/// \file \brief This file implements the boundary condition interface
////////////////////////////////////////////////////////////////////////////////

/// \brief Boundary conditions
namespace boundary {

/// \brief Boundary condition interface
///
/// Single solver BC: apply bc to solver's cells
///
/// Solver/Solver BC: apply bc to active_solver's cells using information from
/// passive_solver
struct Condition {
  void operator()(AnyRange<Ind> cells) const { apply(cells); }
 private:
  virtual void apply(AnyRange<Ind> /*cells*/) const {
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
  /// \brief Constructs a boundary with a custom boundary condition
  template<class Geometry>
  Interface(std::string bcName, std::shared_ptr<Geometry> geometry,
            std::vector<SInd> solverIds,
            std::vector<ConditionHandler> bc)
      : name_(bcName),
        signed_distance([=](const NumA<nd> x){ return geometry->operator()(x); }),
        conditions_(bc),
        solverIds_(solverIds)
  {}

  Interface() = default;
  Interface(const Interface&) = default;
  Interface& operator=(const Interface&) = default;

  std::string name_;
  /// \brief Signed-distance to the boundary
  std::function<Num(const NumA<nd>&)> signed_distance;
  /// \brief Applies boundary condition on boundary cell \p localBCellId of
  /// solver \p s
  std::string name() const { return name_; }
  std::vector<SInd>& solver_ids() { return solverIds_; } // for this solver ids the bc is valid
  const std::vector<SInd>& solver_ids() const { return solverIds_; }
  bool is_valid(const SInd solverId) const {
    return boost::find(solverIds_,solverId) == std::end(solverIds_) ? false : true;
  }
  const Condition& condition(const SInd solverId) const {
    ASSERT(is_valid(solverId), "Invalid BC for solver with it: " << solverId
           << "! Use the is_valid(solverId) function to avoid calling invalid BCs!");
    auto pos = boost::find(solverIds_,solverId) - std::begin(solverIds_);
    return (*conditions_[pos]);
  }
 private:
  std::vector<ConditionHandler> conditions_;
  std::vector<SInd> solverIds_;

};

} // boundary namespace

////////////////////////////////////////////////////////////////////////////////
#endif
