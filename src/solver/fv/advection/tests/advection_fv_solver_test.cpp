/// \brief Test for the FV solver of the Advection equation.
/// Includes:
#include "solver/fv/advection.hpp"
#include "solver/fv/utilities.hpp"
#include "geometry/geometry.hpp"
/// External Includes:
#include "gtest/gtest.h"
/// Options:
#define ENABLE_DBG_ 0
#include "misc/dbg.hpp"
////////////////////////////////////////////////////////////////////////////////

using namespace hom3;

/// General properties:
static const Ind outputInterval = 10;
static const Ind maxNoTimeSteps = 10000;

/// General Physical properties:
static const SInd nd = 2;  ///< #of spatial dimensions

/// Advection Solver properties:

/// Select the physics to solve:
namespace adv_physics = solver::fv::advection;
/// Select the numerical flux:
using flux = adv_physics::flux::local_lax_friedrichs;
/// Select the time integration method:
using time_integration = solver::fv::time_integration::euler_forward;
/// Build a FV solver for advection conduction:
template<SInd nd_>
using AdvSolver = adv_physics::Solver<nd_, flux, time_integration>;
/// Unique identifier for the solver:
static const auto advSolverIdx = SolverIdx{0};

/// Grid properties:

/// Root cell covering the domain [0,1] in each spatial dimension
const auto rootCell = grid::RootCell<nd> {
    NumA<nd>::Constant(0), NumA<nd>::Constant(1)
};
/// Minimum refinement level (for the default minLevel mesh generator):
const SInd minRefLevel = 7;
/// Create a cube grid:
auto test_grid = grid::Grid<nd> {
  grid::helpers::cube::properties<nd>(rootCell, minRefLevel)
};

/// \brief Creates properties for the advection solver:
template<SInd nd>
auto adv_properties(grid::Grid<nd>* grid, const Num timeEnd) {
  using namespace grid::helpers::cube; using namespace io;
  using namespace quantity;
  using InitialDomain = std::function<bool(const NumA<nd>)>;

  const Ind maxNoCells = no_solver_cells_with_gc<nd>(minRefLevel);

  /// \todo Replace with Boost.Log
  std::cerr << "fv max no cells: " << maxNoCells
            << " #leafs: " << no_leaf_nodes<nd>(minRefLevel)
            << " #faces: " << no_faces<nd>()
            << " #nodes/face: "
            << no_nodes_per_face<nd>(minRefLevel) << "\n";

  InitialDomain initialDomain = [](const NumA<nd>) { return true; };

  using Velocity = std::function<NumA<nd>(NumA<nd>)>;

  Velocity radial_velocity = [](const NumA<nd> x) {
    Num u_mag =0.5;
    NumA<nd> x_center = NumA<nd>::Constant(0.5);
    Num radius = (x - x_center).norm();
    const Num angle = 0.5 * math::pi
                      - std::acos((x(0) - x_center(0))/radius);
    NumA<nd> v = NumA<nd>::Zero();
    v(0) = std::abs(radius * u_mag * std::sin(angle + 0.5 * math::pi));
    v(1) = std::abs(radius * u_mag * std::cos(angle + 0.5 * math::pi));
    if (x(1) - x_center(1) > 0.) {
      v(0) *= -1;
    }
    if (x(0) - x_center(0) < 0.) {
      v(1) *= -1;
    }
    return v;
  };

  Properties p;
  insert<grid::Grid<nd>*>  (p, "grid", grid);
  insert<Ind>              (p, "maxNoCells", maxNoCells);
  insert<bool>             (p, "restart", false);
  insert<InitialDomain>    (p, "initialDomain", initialDomain);
  insert<Num>              (p, "CFL", 0.5);

  insert<Num>              (p, "timeEnd", timeEnd);
  insert<Velocity>         (p, "velocity", radial_velocity);
  return p;
}

TEST(advection_fv_solver, one_dimensional_dirichlet) {
  /// Create a solver
  auto advSolver = AdvSolver<nd> {
    advSolverIdx, adv_properties<nd>(&test_grid, 100)
  };

  /// Define an initial condition
  auto constant_ic = [](const NumA<nd> x) {
    Num radius = 0.08;
    NumA<nd> x_center = NumA<nd>::Constant(0.5);
    x_center(0) = 0.25;
    if ((x - x_center).norm() - radius > 0.) {
      return NumA<AdvSolver<nd>::nvars>::Constant(0.1);
    } else {
      return NumA<AdvSolver<nd>::nvars>::Constant(1);
    }
  };
  advSolver.set_initial_condition(constant_ic);

  /// Create Boundary conditions
  auto cBc = adv_physics::bc::Characteristic<AdvSolver<nd>>{advSolver};
  auto bcs = std::make_tuple(cBc, cBc, cBc, cBc, cBc, cBc);
  solver::fv::append_bcs(advSolver, rootCell, bcs);

  /// Run the simulation
  solver::fv::run_solver(test_grid, advSolver, maxNoTimeSteps, outputInterval);
}

////////////////////////////////////////////////////////////////////////////////
#undef ENABLE_DBG_
////////////////////////////////////////////////////////////////////////////////
