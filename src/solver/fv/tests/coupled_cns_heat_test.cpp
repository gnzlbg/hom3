/// \brief Tests for Conjugate Heat-Transfer: Heat-Equation + Navier-Stokes
/// for ideal compressible flow.
////////////////////////////////////////////////////////////////////////////////
/// Includes:
#include "solver/fv/cns.hpp"
#include "solver/fv/cns/coupling_conditions.hpp"
#include "solver/fv/heat.hpp"
#include "solver/fv/heat/coupling_conditions.hpp"
#include "solver/fv/utilities.hpp"
#include "geometry/geometry.hpp"
/// External Includes:
#include "gtest/gtest.h"
/// Options:
#define ENABLE_DBG_ 0
#include "misc/dbg.hpp"
////////////////////////////////////////////////////////////////////////////////
using namespace hom3;

// General properties:
static const Ind outputInterval = 50;
static const Ind maxNoSteps = 40000;

// General Physical properties:
static const SInd nd = 2;
static const Num timeEnd = 10.0;

// CNS Physics:
namespace cns_physics = solver::fv::cns;
/// Select the numerical flux: standard = AUSM + Three-Point stencil
using cns_flux = cns_physics::flux::standard;
/// Select the time integration method:
using cns_time_integration = solver::fv::time_integration::euler_forward;

/// FV CNS Solver:
template<SInd nd_>
using CNSSolver = cns_physics::Solver<nd_, cns_flux, cns_time_integration>;
/// CNS Solver Idx:
static const auto cnsSolverIdx = SolverIdx{0};

/// Boundary/Initial condition: fluid
static const Num jetTRatio = 1.5;

// Heat Physics:
namespace heat_physics = solver::fv::heat;
/// Select the numerical flux:
using heat_flux = heat_physics::flux::three_point;
/// Select the time integration method:
using heat_time_integration = solver::fv::time_integration::euler_forward;

/// FV Heat Solver:
template<SInd nd_>
using HeatSolver = heat_physics::Solver<nd_, heat_flux, heat_time_integration>;
/// Heat Solver Idx:
static const auto heatSolverIdx = SolverIdx{1};

/// Coupling properties:

/// Thermal diffusivity ratio:
static const auto Tr0 = 0.01;

// Grid:
const SInd minRefLevel = 8;

/// Root cell covering the domain [0,1] in each spatial dimension
const auto rootCell = grid::RootCell<nd> {
  NumA<nd>{-7., -7.}, NumA<nd>{7.,7.}
};

/// Create a cube grid:
auto hgrid = grid::Grid<nd> {
  grid::helpers::cube::properties<nd>(rootCell, minRefLevel, 2)
};

/// \brief Makes a cut-off cube
template<SInd nd> auto make_cube
(const NumA<nd> x_center, const NumA<nd> dimensions0, const Num cell_length) {
  const NumA<nd> dimensions1 = dimensions0.array() + (cell_length);
  const NumA<nd> dimensions2 = dimensions0.array() - 0.6 * (cell_length);

  auto cube0
    = geometry::make_geometry<
        geometry::implicit::Square<nd>
      >(x_center, dimensions0);

  auto cube1
    = geometry::make_geometry<
        geometry::implicit::Square<nd>
      >(x_center, dimensions1);

  auto cube2
    = geometry::implicit::adaptors::invert
      (geometry::make_geometry<
         geometry::implicit::Square<nd>
       >(x_center, dimensions2));

  return std::make_tuple(cube0, cube1, cube2);
}

const Num length = 1.0;
const Num h = hgrid.cell_length_at_level(minRefLevel);
const Num x0_cube = -3.0;
auto cubeD = make_cube<nd>(NumA<nd>{x0_cube, 0.0}, NumA<nd>{length, length}, h);

/// \brief Creates properties for the CNS solver
template<SInd nd> auto cns_properties() {
  using namespace grid::helpers::cube; using namespace io;
  using namespace quantity; using namespace unit;
  using InitialDomain = typename CNSSolver<nd>::InitialDomain;

  const Ind maxNoCells = 66308;

  /// \todo Replace with Boost.Log
  std::cerr << "fv max no cells: " << maxNoCells
            << " #leafs: " << no_leaf_nodes<nd>(minRefLevel)
            << " #faces: " << no_faces<nd>()
            << " #nodes/face: "
            << no_nodes_per_face<nd>(minRefLevel) << "\n";

  InitialDomain initialDomain = [&](const NumA<nd> x) {
    return (*std::get<0>(cubeD))(x) > 0. ? true : false;
  };

  Properties p;
  insert<grid::Grid<nd>*>     (p, "grid"         , &hgrid);
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
  insert<Dimensionless>       (p, "Minfty"       , 0.2);

  /// Navier-Stokes properties:
  insert<Dimensionless>       (p, "ReInfty"      , 100);
  insert<Dimensionless>       (p, "Pr0"          , 0.72);
  insert<Temperature>         (p, "Sutherland0"  , 110.4 * kelvin);

  return p;
}

/// \brief Creates properties for the heat solver:
template<SInd nd>
auto heat_properties() {
  using namespace grid::helpers::cube; using namespace io;
  using namespace quantity;
  using InitialDomain = std::function<bool(const NumA<nd>)>;

  const Ind maxNoCells = 396;

  /// \todo Replace with Boost.Log
  std::cerr << "fv max no cells: " << maxNoCells
            << " #leafs: " << no_leaf_nodes<nd>(minRefLevel)
            << " #faces: " << no_faces<nd>()
            << " #nodes/face: "
            << no_nodes_per_face<nd>(minRefLevel) << "\n";

   InitialDomain initialDomain = [&](const NumA<nd> x) {
    return (*std::get<0>(cubeD))(x) < 0. ? true : false;
  };

  Properties p;
  insert<grid::Grid<nd>*>  (p, "grid", &hgrid);
  insert<Ind>              (p, "maxNoCells", maxNoCells);
  insert<bool>             (p, "restart", false);
  insert<InitialDomain>    (p, "initialDomain", initialDomain);
  insert<Num>              (p, "CFL", 0.5);

  insert<Num>                (p, "timeEnd", timeEnd);
  insert<Temperature>        (p, "T_ref", 273.15 * unit::kelvin);
  insert<ThermalDiffusivity> (p, "diffusivity", 1.0 * unit::meter
                                                * unit::meter / unit::second);
  return p;
}

/// \brief Flow past a heat-conducting solid cube
///
///
TEST(conjugate_heat_transfer_fv_solver, flow_past_solid_cube) {
  using namespace grid::helpers::cube;

  /// Create CNS solver:
  auto cnsSolver = CNSSolver<nd> {
    cnsSolverIdx, cns_properties<nd>()
  };

  /// CNS initial condition:
  auto ic = [&](const NumA<nd>) {
    using V = CNSSolver<nd>::V;
    NumA<cnsSolver.nvars> pv = NumA<cnsSolver.nvars>::Zero();
    pv(V::u(0))  = cnsSolver.quantities.u_infinity();
    pv(V::p())   = cnsSolver.quantities.p_infinity();
    pv(V::rho()) = cnsSolver.quantities.rho_infinity();
    return cnsSolver.cv(pv);
  };
  cnsSolver.set_initial_condition(ic);

  /// CNS far-field boundary conditions:
  auto nBc = cns_physics::bc::Neumann<CNSSolver<nd>>{cnsSolver};

  /// CNS inflow-outflow boundary conditions:
  auto free_stream = [&](const CellIdx, const SInd v) {
        using V = CNSSolver<nd>::V;
    NumA<cnsSolver.nvars> pv = NumA<cnsSolver.nvars>::Zero();
    pv(V::u(0))  = cnsSolver.quantities.u_infinity();
    pv(V::p())   = cnsSolver.quantities.p_infinity();
    pv(V::rho()) = cnsSolver.quantities.rho_infinity();
    return pv(v);
  };

  auto iBc = cns_physics::bc::external::Inflow<CNSSolver<nd>> {
    cnsSolver, free_stream
  };
  auto oBc = cns_physics::bc::external::Outflow<CNSSolver<nd>> {
    cnsSolver, free_stream
  };

  /// CNS freestream boundaries:
  auto bcs = std::make_tuple(iBc, oBc, nBc, nBc);
  solver::fv::append_bcs(cnsSolver, hgrid.root_cell(), bcs);

  /// Create Heat solver:
  auto heatSolver = HeatSolver<nd> {
    heatSolverIdx, heat_properties<nd>()
  };

  /// Heat initial condition:
  auto constant_ic = [&](const NumA<nd>) {
    return NumA<heatSolver.nvars>::Constant(cnsSolver.quantities.T_infinity());
  };
  heatSolver.set_initial_condition(constant_ic);

  /// Heat coupling condition:
  auto hwBc = heat_physics::bc::coupling::CNS<HeatSolver<nd>, CNSSolver<nd>> {
     heatSolver, cnsSolver
  };
  heatSolver.append_bc(solver::fv::bc::Interface<nd> {
      "cubeHeat", std::get<2>(cubeD), heatSolver, hwBc
  });

  /// CNS coupling condition:
  auto cwBc = cns_physics::bc::coupling::Heat<CNSSolver<nd>, HeatSolver<nd>> {
    cnsSolver, heatSolver
  };
  cnsSolver.append_bc(solver::fv::bc::Interface<nd> {
      "cubeCNS", std::get<1>(cubeD), cnsSolver, cwBc
  });

  /// Initialize:
  initialize(hgrid, heatSolver, cnsSolver);
  write_domains(hgrid, heatSolver, cnsSolver);

  /// Main loop:
  while (!solver_finished(maxNoSteps, cnsSolver, heatSolver)) {
    solver::fv::coupling::apply(Tr0, cnsSolver, heatSolver);
    write_timestep(cnsSolver, heatSolver);

    cnsSolver.solve();
    heatSolver.solve();

    write_output(outputInterval, cnsSolver, heatSolver);

    if(solution_diverged(cnsSolver, heatSolver)) {
      write_domains(cnsSolver, heatSolver);
      TERMINATE("Solution diverged!");
    }
  }
  write_domains(cnsSolver, heatSolver);

}

////////////////////////////////////////////////////////////////////////////////
#undef ENABLE_DBG_
////////////////////////////////////////////////////////////////////////////////
