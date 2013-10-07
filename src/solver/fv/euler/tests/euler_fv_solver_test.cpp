/// \brief Tests for the FV solver of the Euler equations.
////////////////////////////////////////////////////////////////////////////////
/// Includes:
#include "solver/fv/euler.hpp"
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
static const Ind outputInterval = 50;
static const Ind maxNoTimeSteps = 10000;

/// General Physical properties:
static const SInd nd = 3;

/// Euler Solver properties:
static const auto eulerSolverIdx = SolverIdx{0};

/// Euler FV Solver:
namespace euler_physics = solver::fv::euler;

using flux = euler_physics::flux::ausm;
using time_integration = solver::fv::time_integration::euler_forward;

template<SInd nd_>
using EulerSolver = euler_physics::Solver<nd_, flux, time_integration>;

using V = typename EulerSolver<nd>::V;

/// Grid:
const auto rootCell = grid::RootCell<nd> {
  NumA<nd>::Constant(0), NumA<nd>::Constant(1)
};
const SInd minRefLevel = 5;
auto test_grid = grid::Grid<nd> {
  grid::helpers::cube::properties<nd>(rootCell, minRefLevel)
};

/// \brief Creates properties for the Euler solver
template<SInd nd>
auto euler_properties(grid::Grid<nd>* grid, const Num timeEnd) {
  using namespace grid::helpers::cube; using namespace io;
  using namespace quantity; using namespace unit;
  using InitialDomain = std::function<bool(const NumA<nd>)>;

  const Ind maxNoCells = no_solver_cells_with_gc<nd>(minRefLevel);

  /// \todo Replace with Boost.Log
  std::cerr << "fv max no cells: " << maxNoCells
            << " #leafs: " << no_leaf_nodes<nd>(minRefLevel)
            << " #faces: " << no_faces<nd>()
            << " #nodes/face: "
            << no_nodes_per_face<nd>(minRefLevel) << "\n";

  InitialDomain initialDomain = [](const NumA<nd>) { return true; };

  Properties p;
  insert<grid::Grid<nd>*>     (p, "grid"         , grid);
  insert<Ind>                 (p, "maxNoCells"   , maxNoCells);
  insert<bool>                (p, "restart"      , false);
  insert<InitialDomain>       (p, "initialDomain", initialDomain);
  insert<Num>                 (p, "timeEnd"      , timeEnd);
  insert<Dimensionless>       (p, "gamma"        , 1.4);
  insert<Dimensionless>       (p, "Minfty"       , 0.1);
  insert<Num>                 (p, "CFL"          , 0.2);
  insert<SpecificGasConstant> (p, "Rspecific"    , 287.058 * joule
                                                   / (kilogram * kelvin));
  insert<Temperature>         (p, "T0"           , 273.15 * kelvin);
  insert<Density>             (p, "rho0"         , 1.225 * kilogram
                                                        / cubic_meter);
  return p;
}


/// \test Constant initial condition with Neumann boundary conditions
///
/// Note: The solution is the initial condition for all $t$
TEST(euler_fv_solver, constant_ic) {
  using namespace grid::helpers::cube;

  /// Create solver
  auto eulerSolver = EulerSolver<nd> {
    eulerSolverIdx, euler_properties<nd>(&test_grid, 0.25)
  };

  /// Define an initial condition
  auto constant_ic = [&](const NumA<nd>) {
    NumA<eulerSolver.nvars> pvars = NumA<eulerSolver.nvars>::Zero();
    pvars(V::rho()) = 1.0;
    for (SInd d = 0; d < nd; ++d) {
      pvars(V::u(d)) = 0.5;
    }
    pvars(V::p()) = 1.0;
    return eulerSolver.cv(pvars);
  };
  eulerSolver.set_initial_condition(constant_ic);

  /// Create boundary conditions
  auto nBc = euler_physics::bc::Neumann<EulerSolver<nd>>(eulerSolver);
  solver::fv::append_bcs(eulerSolver, test_grid.root_cell(),
                         make_conditions<nd>(nBc));

  solver::fv::run_solver(test_grid, eulerSolver, maxNoTimeSteps,
                         outputInterval);
}

/// \test Modified Sod's test
///
/// This is a really easy problem. If something goes wrong here,
/// something fundamental is wrong.
///
///     rho  |   u   |  p
/// L:  1.0  | 0.75  | 1.0
/// R: 0.125 | 0.0   | 0.1
///
/// see Riemann Solvers and Numerical Methods for Fluid Dynamics 3rd
/// Edition p. 225
///
TEST(euler_fv_solver, sod_shock_tube_ic) {
  using namespace grid::helpers::cube;

  /// Create solver
  auto eulerSolver = EulerSolver<nd> {
    eulerSolverIdx, euler_properties<nd>(&test_grid, 0.2)
  };

  /// Define initial condition
  const SInd dir = 0;
  const Num angle = 0;
  const Num x0 = 0.3;
  const Num rhoL = 1.0;
  const Num umagL = 0.75;
  const Num pL = 1.0;
  const Num rhoR = 0.125;
  const Num umagR = 0.0;
  const Num pR = 0.1;
  auto sod_shock_tube_ic
    = euler_physics::ic::shock_tube<nd>
      (dir, angle, x0, rhoL, umagL, pL, rhoR, umagR, pR);
  eulerSolver.set_initial_condition(sod_shock_tube_ic);

  /// Create boundary conditions
  auto nBc = euler_physics::bc::Neumann<EulerSolver<nd>>(eulerSolver);
  solver::fv::append_bcs(eulerSolver, test_grid.root_cell(),
                         make_conditions<nd>(nBc));

  solver::fv::run_solver(test_grid, eulerSolver, maxNoTimeSteps,
                         outputInterval);
}

/// \test 123 problem
///
/// This is a tough problem that asses the performance for low-density flows.
/// The checks for negative density/pressure should fail before anything goes
/// wrong here.
///
///     rho  |   u   |  p
/// L:  1.0  | -2.0  | 0.4
/// R:  1.0  |  2.0  | 0.4
///
/// see Riemann Solvers and Numerical Methods for Fluid Dynamics 3rd
/// Edition p. 225
///
TEST(euler_fv_solver, problem123_ic) {
  using namespace grid::helpers::cube;

  /// Create solver
  auto eulerSolver = EulerSolver<nd> {
    eulerSolverIdx, euler_properties<nd>(&test_grid, 0.15)
  };

  /// Define initial condition
  const SInd dir = 0;
  const Num angle = 0;
  const Num x0 = 0.5;
  const Num rhoL = 1.0;
  const Num umagL = -2.0;
  const Num pL = 0.4;
  const Num rhoR = 1.0;
  const Num umagR = 2.0;
  const Num pR = 0.4;
  auto problem_123_ic
    = euler_physics::ic::shock_tube<nd>
      (dir, angle, x0, rhoL, umagL, pL, rhoR, umagR, pR);
  eulerSolver.set_initial_condition(problem_123_ic);

  /// Create boundary conditions
  auto nBc = euler_physics::bc::Neumann<EulerSolver<nd>>(eulerSolver);
  solver::fv::append_bcs(eulerSolver, test_grid.root_cell(),
                         make_conditions<nd>(nBc));

  solver::fv::run_solver(test_grid, eulerSolver, maxNoTimeSteps,
                         outputInterval);
}

/// \test Explosion
///
/// Discontinuity $\mathbf{x}_0$ is a sphere centered
/// at $\mathbf{x}_c = (0.5, 0.5, 0.5)^T of radius $r = 0.2$.
///
///           rho   |  u  |  p
/// Inside:   1.0   | 0.0 | 1.0
/// Outside:  0.125 | 0.0 | 0.1
///
TEST(euler_fv_solver, explosion_ic) {
  using namespace grid::helpers::cube;

  /// Create solver
  auto eulerSolver = EulerSolver<nd> {
    eulerSolverIdx, euler_properties<nd>(&test_grid, 0.25)
  };

  /// Define initial condition
  const Num radius = 0.2;
  auto explosion_ic = [&](const NumA<nd> x) {
    const Num R = radius;
    NumA<nd> xE;
    xE(0) = 0.5; xE(1) = 0.5;
    if (nd == 3) {
      xE(2) = 0.5;
    }

    NumA<eulerSolver.nvars> pvars;
    if ((x - xE).norm() < R) {  // inside
      pvars(V::rho()) = 1.0;
      for (SInd d = 0; d < nd; ++d) {
        pvars(V::u(d)) = 0;
      }
      pvars(V::p()) = 1.0;
      } else {  // outside
      pvars(V::rho()) = 0.125;
      for (SInd d = 0; d < nd; ++d) {
        pvars(V::u(d)) = 0;
      }
      pvars(V::p()) = 0.1;
    }
    return eulerSolver.cv(pvars);
  };
  eulerSolver.set_initial_condition(explosion_ic);

  /// Create boundary conditions
  auto nBc = euler_physics::bc::Neumann<EulerSolver<nd>>(eulerSolver);
  solver::fv::append_bcs(eulerSolver, test_grid.root_cell(),
                         make_conditions<nd>(nBc));

  solver::fv::run_solver(test_grid, eulerSolver, maxNoTimeSteps,
                         outputInterval);
}

/// \test Isentropic Vortex (Euler-equations 2D)
TEST(euler_fv_solver, isentropic_vortex_ic) {
  using namespace grid::helpers::cube;

  /// Create grid
  static const SInd nd = 2;
  const auto rootCell_2d = grid::RootCell<nd> {
    NumA<nd>::Constant(-8.0), NumA<nd>::Constant(8.0)
  };
  auto test_grid_2d = grid::Grid<nd>{properties<nd>(rootCell_2d, minRefLevel)};

  /// Create solver
  auto eulerSolver = EulerSolver<nd> {
    eulerSolverIdx, euler_properties<nd>(&test_grid_2d, 2.0)
  };

  /// Define initial condition
  auto isentropic_vortex_ic = [](const NumA<nd> x) {
    return euler_physics::ic::isentropic_vortex<nd>(x, 0);
  };
  eulerSolver.set_initial_condition(isentropic_vortex_ic);

  /// Create boundary conditions
  auto nBc = euler_physics::bc::IsentropicVortex<EulerSolver<nd>>(eulerSolver);
  solver::fv::append_bcs(eulerSolver, test_grid_2d.root_cell(),
                         make_conditions<nd>(nBc));

  solver::fv::run_solver(test_grid_2d, eulerSolver, maxNoTimeSteps,
                         outputInterval);
}

////////////////////////////////////////////////////////////////////////////////
#undef ENABLE_DBG_
////////////////////////////////////////////////////////////////////////////////
