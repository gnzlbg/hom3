/// \brief Tests for Euler-Euler coupling.
////////////////////////////////////////////////////////////////////////////////
/// Includes:
#include "solver/fv/euler.hpp"
#include "solver/fv/euler/coupling_conditions.hpp"
#include "solver/fv/coupling.hpp"
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
static const Ind maxNoSteps = 40000;

/// General Physical properties:
static const SInd nd = 3;
static const Num timeEnd = 0.25;

/// Euler Physics:
namespace euler_physics = solver::fv::euler;
/// Select the numerical flux: AUSM
using flux = euler_physics::flux::ausm;
/// Select the time integration method:
using time_integration = solver::fv::time_integration::euler_forward;

/// FV Euler Solver:
template<SInd nd_>
using EulerSolver = euler_physics::Solver<nd_, flux, time_integration>;
// Euler solver indices:
static const auto eulerSolverIdx0 = SolverIdx{0};
static const auto eulerSolverIdx1 = SolverIdx{1};

/// Grid:
const SInd minRefLevel = 5;

/// Root cell covering the domain [0,1] in each spatial dimension
const auto rootCell = grid::RootCell<nd> {
  NumA<nd>::Constant(0), NumA<nd>::Constant(1)
};

/// Create a cube grid:
auto test_grid = grid::Grid<nd> {
  grid::helpers::cube::properties<nd>(rootCell, minRefLevel, 2)
};

const Num length = 0.2;
const Num h = test_grid.cell_length_at_level(minRefLevel);
auto cubeD = geometry::make_cube<nd>
             (NumA<nd>::Constant(0.5), NumA<nd>::Constant(length), h, 0.5);

/// \brief Creates properties for the Euler solvers
template<SInd nd, class InitD> auto euler_properties(InitD&& id) {
  using namespace grid::helpers::cube; using namespace io;
  using namespace quantity; using namespace unit;
  using InitialDomain = typename EulerSolver<nd>::InitialDomain;

  const Ind maxNoCells = no_solver_cells_with_gc<nd>(minRefLevel) * 2;

  /// \todo Replace with Boost.Log
  std::cerr << "fv max no cells: " << maxNoCells
            << " #leafs: " << no_leaf_nodes<nd>(minRefLevel)
            << " #faces: " << no_faces<nd>()
            << " #nodes/face: "
            << no_nodes_per_face<nd>(minRefLevel) << "\n";

  InitialDomain initialDomain = id;

  Properties p;
  insert<grid::Grid<nd>*>     (p, "grid"         , &test_grid);
  insert<Ind>                 (p, "maxNoCells"   , maxNoCells);
  insert<bool>                (p, "restart"      , false);
  insert<InitialDomain>       (p, "initialDomain", initialDomain);
  insert<Num>                 (p, "timeEnd"      , timeEnd);
  insert<Num>                 (p, "CFL"          , 0.2);

  /// Euler properties:
  insert<Temperature>         (p, "T0"           , 273.15 * kelvin);
  insert<Density>             (p, "rho0"         , 1.225 * kilogram
                                                        / cubic_meter);
  insert<SpecificGasConstant> (p, "Rspecific"    , 287.058 * joule
                                                   / (kilogram * kelvin));
  insert<Dimensionless>       (p, "gamma"        , 1.4);
  insert<Dimensionless>       (p, "Minfty"       , 0.2);

  return p;
}

TEST(euler_fv_solver_coupling, euler_euler_coupling) {
  using namespace grid::helpers::cube;

  /// Create Euler solvers:
  auto eulerSolver0 = EulerSolver<nd> {
    eulerSolverIdx0, euler_properties<nd>([&](const NumA<nd> x) {
    return (*std::get<0>(cubeD))(x) > 0. ? true : false;
      })
  };
  auto eulerSolver1 = EulerSolver<nd> {
    eulerSolverIdx1,
    euler_properties<nd>([&](const NumA<nd> x) {
    return (*std::get<0>(cubeD))(x) < 0. ? true : false;
      })
  };

  /// Set initial condition:
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
  eulerSolver0.set_initial_condition(sod_shock_tube_ic);
  eulerSolver1.set_initial_condition(sod_shock_tube_ic);

  /// Transmissive boundary conditions:
  auto nBc = euler_physics::bc::Neumann<EulerSolver<nd>>(eulerSolver0);
  solver::fv::append_bcs(eulerSolver0, test_grid.root_cell(),
                         make_conditions<nd>(nBc));

  /// Coupling condition:
  auto cBc0 = euler_physics::bc::coupling::Euler<EulerSolver<nd>> {
    eulerSolver0, eulerSolver1
  };
  eulerSolver0.append_bc(solver::fv::bc::Interface<nd> {
      "cube0", std::get<1>(cubeD), eulerSolver0, cBc0
  });
  auto cBc1 = euler_physics::bc::coupling::Euler<EulerSolver<nd>> {
    eulerSolver1, eulerSolver0
  };
  eulerSolver1.append_bc(solver::fv::bc::Interface<nd> {
      "cube1", std::get<2>(cubeD), eulerSolver1, cBc1
  });

  /// Initialize:
  initialize(test_grid, eulerSolver0, eulerSolver1);
  write_domains(test_grid, eulerSolver0, eulerSolver1);

  /// Main loop:
  while (!solver_finished(maxNoSteps, eulerSolver0, eulerSolver1)) {
    solver::fv::coupling::apply(eulerSolver0, eulerSolver1);
    write_timestep(eulerSolver0, eulerSolver1);

    eulerSolver0.solve();
    eulerSolver1.solve();

    write_output(outputInterval, eulerSolver0, eulerSolver1);

    if (solution_diverged(eulerSolver0, eulerSolver1)) {
      write_domains(eulerSolver0, eulerSolver1);
      TERMINATE("Solution diverged!");
    }
  }
  write_domains(eulerSolver0, eulerSolver1);
}

////////////////////////////////////////////////////////////////////////////////
#undef ENABLE_DBG_
////////////////////////////////////////////////////////////////////////////////
