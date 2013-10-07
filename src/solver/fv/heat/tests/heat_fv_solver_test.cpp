/// \brief Test for the FV solver of the Heat equation.
/// Includes:
#include "solver/fv/heat.hpp"
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
static const Ind maxNoTimeSteps = 1000;

/// General Physical properties:
static const SInd nd = 3;  ///< #of spatial dimensions

/// Heat Solver properties:

/// Select the physics to solve:
namespace heat_physics = solver::fv::heat;
/// Select the numerical flux:
using flux = heat_physics::flux::three_point;
/// Select the time integration method:
using time_integration = solver::fv::time_integration::euler_forward;
/// Build a FV solver for heat conduction:
template<SInd nd_>
using HeatSolver = heat_physics::Solver<nd_, flux, time_integration>;
/// Unique identifier for the solver:
static const auto heatSolverIdx = SolverIdx{0};

/// Grid properties:

/// Root cell covering the domain [0,1] in each spatial dimension
const auto rootCell = grid::RootCell<nd> {
    NumA<nd>::Constant(0), NumA<nd>::Constant(1)
};
/// Minimum refinement level (for the default minLevel mesh generator):
const SInd minRefLevel = 4;
/// Create a cube grid:
auto test_grid = grid::Grid<nd> {
  grid::helpers::cube::properties<nd>(rootCell, minRefLevel)
};

/// \brief Creates properties for the heat solver:
template<SInd nd>
auto heat_properties(grid::Grid<nd>* grid, const Num timeEnd) {
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

  Properties p;
  insert<grid::Grid<nd>*>  (p, "grid", grid);
  insert<Ind>              (p, "maxNoCells", maxNoCells);
  insert<bool>             (p, "restart", false);
  insert<InitialDomain>    (p, "initialDomain", initialDomain);
  insert<Num>              (p, "CFL", 1.0);

  insert<Num>                (p, "timeEnd", timeEnd);
  insert<Temperature>        (p, "T_ref", 273.15 * unit::kelvin);
  insert<ThermalDiffusivity> (p, "diffusivity", (1.0 * unit::meter * unit::meter
                                                      / unit::second));
  return p;
}

TEST(heat_fv_solver, one_dimensional_dirichlet) {
  /// Create a solver
  auto heatSolver =   HeatSolver<nd> {
    heatSolverIdx, heat_properties<nd>(&test_grid, 0.25)
  };

  /// Define an initial condition
  auto constant_ic = [](const NumA<nd>) {
    return NumA<HeatSolver<nd>::nvars>::Zero();
  };
  heatSolver.set_initial_condition(constant_ic);

  /// Create Boundary conditions
  auto dBc = heat_physics::bc::Dirichlet<HeatSolver<nd>>{heatSolver, 1.0};
  auto nBc = heat_physics::bc::Neumann<HeatSolver<nd>>{heatSolver, 0.0};
  auto bcs = std::make_tuple(dBc, dBc, nBc, nBc, nBc, nBc);
  solver::fv::append_bcs(heatSolver, rootCell, bcs);

  /// Run the simulation
  solver::fv::run_solver(test_grid, heatSolver, maxNoTimeSteps, outputInterval);
}

////////////////////////////////////////////////////////////////////////////////
#undef ENABLE_DBG_
////////////////////////////////////////////////////////////////////////////////
