/// \brief Tests for the FV solver of the Compressible Navier-Stokes equations.
////////////////////////////////////////////////////////////////////////////////
/// Includes:
#include "solver/fv/cns.hpp"
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
static const Ind outputInterval = 20;
static const Ind maxNoTimeSteps = 40000;

/// General Physical properties:
static const SInd nd = 2;  ///< #of spatial dimensions

/// CNS Solver properties:

/// Select the physics to solve:
namespace cns_physics = solver::fv::cns;
/// Select the numerical flux: standard = AUSM + Three-Point stencil
using flux = cns_physics::flux::standard;
/// Select the time integration method:
using time_integration = solver::fv::time_integration::euler_forward;
/// Build a FV solver for the CNS equations:
template<SInd nd_>
using CNSSolver = cns_physics::Solver<nd_, flux, time_integration>;
/// Unique identifier for the solver:
static const auto cnsSolverIdx = SolverIdx{0};
/// Alias the access to the solver variables:
using V = typename CNSSolver<nd>::V;
/// #of variables
static const auto nvars = CNSSolver<nd>::nvars;

/// Grid properties:

/// Minimum refinement level (for the default minLevel mesh generator):
const SInd minRefLevel = 8;

/// \brief Creates properties for the Euler solver
template<SInd nd, class InitD>
auto cns_properties
(grid::Grid<nd>* grid, const Num timeEnd, const Num ReInfinity,
  InitD initialDomain = [](const NumA<nd>) { return true; },
  const Num MInfinity = 0.1) {
  using namespace grid::helpers::cube; using namespace io;
  using namespace quantity; using namespace unit;
  using InitialDomain = typename CNSSolver<nd>::InitialDomain;

  const Ind maxNoCells = 1.2 * no_solver_cells_with_gc<nd>(minRefLevel);

  /// \todo Replace with Boost.Log
  std::cerr << "fv max no cells: " << maxNoCells
            << " #leafs: " << no_leaf_nodes<nd>(minRefLevel)
            << " #faces: " << no_faces<nd>()
            << " #nodes/face: "
            << no_nodes_per_face<nd>(minRefLevel) << "\n";

  Properties p;
  insert<grid::Grid<nd>*>     (p, "grid"         , grid);
  insert<Ind>                 (p, "maxNoCells"   , maxNoCells);
  insert<bool>                (p, "restart"      , false);
  insert<InitialDomain>       (p, "initialDomain", initialDomain);
  insert<Num>                 (p, "timeEnd"      , timeEnd);
  insert<Num>                 (p, "CFL"          , 0.4);
  insert<Num>                 (p, "CFL_viscous"  , 0.4);

  /// Euler properties:
  insert<Temperature>         (p, "T0"           , 273.15 * kelvin);
  insert<Density>             (p, "rho0"         , 1.225 * kilogram
                                                        / cubic_meter);
  insert<SpecificGasConstant> (p, "Rspecific"    , 287.058 * joule
                                                   / (kilogram * kelvin));
  insert<Dimensionless>       (p, "gamma"        , 1.4);
  insert<Dimensionless>       (p, "Minfty"       , MInfinity);

  /// Navier-Stokes properties:
  insert<Dimensionless>       (p, "ReInfty"      , ReInfinity);
  insert<Dimensionless>       (p, "Pr0"          , 0.72);
  insert<Temperature>         (p, "Sutherland0"  , 110.4 * kelvin);

  return p;
}

TEST(cns_fv_solver, constant_ic) {
  using namespace grid::helpers::cube;

  /// Root cell covering the domain [0,1] in each spatial dimension
  const auto rootCell = grid::RootCell<nd> {
    NumA<nd>{0., 0.}, NumA<nd>{1., 1.}
  };

  /// Create a cube grid:
  auto grid = grid::Grid<nd> {
    grid::helpers::cube::properties<nd>(rootCell, minRefLevel)
  };

  /// Create solver
  auto cnsSolver = CNSSolver<nd> {
    cnsSolverIdx,
    cns_properties<nd>(&grid, 0.25, 100., [&](const NumA<nd>) {
        return true;
      })
  };

  /// Define an initial condition
  auto constant_ic = [&](const NumA<nd>) {
    NumA<CNSSolver<nd>::nvars> pvars = NumA<CNSSolver<nd>::nvars>::Zero();
    pvars(V::rho()) = cnsSolver.quantities.rho_infinity();
    pvars(V::u(0)) = cnsSolver.quantities.u_infinity();
    pvars(V::p()) = cnsSolver.quantities.p_infinity();
    return cnsSolver.cv(pvars);
  };
  cnsSolver.set_initial_condition(constant_ic);

  /// Create boundary conditions
  auto nBc = cns_physics::bc::Neumann<CNSSolver<nd>>(cnsSolver);
  solver::fv::append_bcs(cnsSolver, grid.root_cell(),
                         make_conditions<nd>(nBc));

  solver::fv::run_solver(grid, cnsSolver, maxNoTimeSteps,
                         outputInterval);
}

/// \test Modified Sod's test
///
/// Mild test: an error here means something fundamental is wrong.
///
///     rho  |   u   |  p
/// L:  1.0  | 0.75  | 1.0
/// R: 0.125 | 0.0   | 0.1
///
/// see Riemann Solvers and Numerical Methods for Fluid Dynamics 3rd
/// Edition p. 225
///
TEST(cns_fv_solver, sod_shock_tube_ic) {
  using namespace grid::helpers::cube;

  /// Root cell covering the domain [0,1] in each spatial dimension
  const auto rootCell = grid::RootCell<nd> {
    NumA<nd>{0., 0.}, NumA<nd>{1., 1.}
  };

  /// Create a cube grid:
  auto grid = grid::Grid<nd> {
    grid::helpers::cube::properties<nd>(rootCell, minRefLevel)
  };

  /// Create solver
  auto cnsSolver = CNSSolver<nd> {
    cnsSolverIdx,
    cns_properties<nd>(&grid, 0.25, 1000., [&](const NumA<nd>) {
        return true;
      })
  };

  /// Define initial condition
  const auto rho_infinity = cnsSolver.quantities.rho_infinity();
  const auto u_infinity   = cnsSolver.quantities.u_infinity();
  const auto p_infinity   = cnsSolver.quantities.p_infinity();
  const SInd dir  = 0;
  const Num angle = 0;
  const Num x0    = 0.3;
  const Num rhoL  = 1.0   * rho_infinity;
  const Num umagL = 0.75  * u_infinity;
  const Num pL    = 1.0   * p_infinity;
  const Num rhoR  = 0.125 * rho_infinity;
  const Num umagR = 0.0   * u_infinity;
  const Num pR    = 0.1   * p_infinity;
  auto sod_shock_tube_ic
    = solver::fv::euler::ic::shock_tube<nd>
      (dir, angle, x0, rhoL, umagL, pL, rhoR, umagR, pR);
  cnsSolver.set_initial_condition(sod_shock_tube_ic);

  /// Create boundary conditions
  auto nBc = cns_physics::bc::Neumann<CNSSolver<nd>>(cnsSolver);
  solver::fv::append_bcs(cnsSolver, grid.root_cell(),
                         make_conditions<nd>(nBc));

  solver::fv::run_solver(grid, cnsSolver, maxNoTimeSteps,
                         outputInterval);
}

template<SInd nd> auto make_cube
(const NumA<nd> x_center, const NumA<nd> dimensions, const Num cell_length) {
  const NumA<nd> dimensions2 = dimensions.array() + (cell_length);
  std::cerr << "cell length: " << cell_length << "\n";
  auto cube0
    = geometry::make_geometry<
        geometry::implicit::Square<nd>
      >(x_center, dimensions);

  auto cube1
    = geometry::make_geometry<
        geometry::implicit::Square<nd>
      >(x_center, dimensions2);

  return std::make_tuple(cube0, cube1);
}

/// \brief Flow past a cube testcase
///
/// Reynolds number: Re = 200
/// Initial condition: free-stream values
/// BCs: free-stream inflow/outflow + Neumann walls + AdiabaticNoSlip
TEST(cns_fv_solver, flow_past_cube) {
  using namespace grid::helpers::cube;

  /// Root cell covering the domain [0,1] in each spatial dimension
  const auto rootCell = grid::RootCell<nd> {
    NumA<nd>{-7., -7.}, NumA<nd>{7., 7.}
  };

  /// Create a cube grid:
  auto grid = grid::Grid<nd> {
    grid::helpers::cube::properties<nd>(rootCell, minRefLevel)
  };

  const Num length = 1.0;
  const Num x0_cube = -3.0;
  auto cube = make_cube<nd>(NumA<nd>{x0_cube, 0.0}, NumA<nd>{length, length},
                              grid.cell_length_at_level(minRefLevel));

  /// Create solver
  auto cnsSolver = CNSSolver<nd> {
    cnsSolverIdx,
    cns_properties<nd>(&grid, 100.0, 200., [&](const NumA<nd> x) {
        return (*std::get<0>(cube))(x) > 0. ? true : false;
      })
  };

  /// Define initial condition
  auto ic = [&](const NumA<nd> x) {
    NumA<nvars> pv = NumA<nvars>::Zero();
    pv(V::u(0))  = cnsSolver.quantities.u_infinity();
    pv(V::p())   = cnsSolver.quantities.p_infinity();
    pv(V::rho()) = cnsSolver.quantities.rho_infinity();

    if (x(0) < x0_cube - length && x(0) > x0_cube - 2. * length
       && x(1) > 0. && x(1) < length) {
      pv(V::rho()) = pv(V::rho()) * 0.9;  // density bump
    }

    return cnsSolver.cv(pv);
  };
  cnsSolver.set_initial_condition(ic);

  /// Create far-field boundary conditions
  auto nBc = cns_physics::bc::Neumann<CNSSolver<nd>>(cnsSolver);

  /// Create inflow-outflow boundary conditions
  auto free_stream = [&](const CellIdx, const SInd v) {
    NumA<nvars> pv = NumA<nvars>::Zero();
    pv(V::u(0))  = cnsSolver.quantities.u_infinity();
    pv(V::p())   = cnsSolver.quantities.p_infinity();
    pv(V::rho()) = cnsSolver.quantities.rho_infinity();
    return pv(v);
  };
  auto iBc = cns_physics::bc::external::Inflow<CNSSolver<nd>>
             (cnsSolver, free_stream);
  auto oBc = cns_physics::bc::external::Outflow<CNSSolver<nd>>
             (cnsSolver, free_stream);

  /// Free-stream boundaries:
  auto bcs = std::make_tuple(iBc, oBc, nBc, nBc);
  solver::fv::append_bcs(cnsSolver, grid.root_cell(), bcs);

  /// Create cube solid wall condition:
  auto wBc = cns_physics::bc::wall::AdiabaticNoSlip<CNSSolver<nd>>(cnsSolver);
  cnsSolver.append_bc(solver::fv::bc::Interface<nd> {
      "cube", std::get<1>(cube), cnsSolver, wBc
  });

  solver::fv::run_solver(grid, cnsSolver, maxNoTimeSteps,
                         outputInterval);
}

////////////////////////////////////////////////////////////////////////////////
#undef ENABLE_DBG_
////////////////////////////////////////////////////////////////////////////////
