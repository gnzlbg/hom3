/// \brief Tests for the Coupled MB FV Conjugate Heat Transfer solver.
////////////////////////////////////////////////////////////////////////////////
/// Includes:
#include "solver/fv/cns.hpp"
#include "solver/fv/cns/moving_boundary_conditions.hpp"
#include "solver/fv/cns/coupling_conditions.hpp"
#include "solver/fv/heat.hpp"
#include "solver/fv/heat/moving_boundary_conditions.hpp"
#include "solver/fv/heat/coupling_conditions.hpp"
#include "solver/fv/utilities.hpp"
#include "geometry/geometry.hpp"
/// External Includes:
#include "misc/test.hpp"
/// Options:
#define ENABLE_DBG_ 0
#include "misc/dbg.hpp"
////////////////////////////////////////////////////////////////////////////////

using namespace hom3;

/// General properties:
static const Ind outputInterval = 25;
static const Ind maxNoTimeSteps = 400000;
static const Ind timeEnd = 3000;

/// General Physical properties:
static const SInd nd = 2;  ///< #of spatial dimensions
static const Num T_scaling = 2.0;

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

/// Alias the access to the solver variables:
using CNSV = typename CNSSolver<nd>::V;
/// #of variables
static const auto cns_nvars = CNSSolver<nd>::nvars;

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

/// Alias the access to the solver variables:
using HV = typename HeatSolver<nd>::V;
/// #of variables
static const auto heat_nvars = HeatSolver<nd>::nvars;

/// Coupling properties:

/// Thermal diffusivity ratio:
static const Num Tr0 = 0.0001;

/// Grid properties:

/// Minimum refinement level (for the default minLevel mesh generator):
const SInd minRefLevel = 8;

/// Root cell covering the domain [0,1] in each spatial dimension
const auto rootCell = grid::RootCell<nd> {
  NumA<nd>{-7., -7.}, NumA<nd>{7., 7.}
};

/// Create a cube grid:
auto hgrid = grid::Grid<nd> {
  grid::helpers::cube::properties<nd>(rootCell, minRefLevel, 2)
};


/// \brief Creates properties for the Euler solver
template<SInd nd>
auto cns_properties() {
  using namespace grid::helpers::cube; using namespace io;
  using namespace quantity; using namespace unit;

  const Ind maxNoCells = 1.2 * no_solver_cells_with_gc<nd>(minRefLevel);

  /// \todo Replace with Boost.Log
  std::cerr << "fv max no cells: " << maxNoCells
            << " #leafs: " << no_leaf_nodes<nd>(minRefLevel)
            << " #faces: " << no_faces<nd>()
            << " #nodes/face: "
            << no_nodes_per_face<nd>(minRefLevel) << "\n";

using InitialDomain = typename CNSSolver<nd>::InitialDomain;
  InitialDomain id = [&](NumA<nd>) { return true; };

  Properties p;
  insert<grid::Grid<nd>*>     (p, "grid"         , &hgrid);
  insert<Ind>                 (p, "maxNoCells"   , maxNoCells);
  insert<bool>                (p, "restart"      , false);
  insert<Num>                 (p, "timeEnd"      , timeEnd);
  insert<Num>                 (p, "CFL"          , 0.5);
  insert<Num>                 (p, "CFL_viscous"  , 0.4);

  insert<InitialDomain>       (p, "initialDomain", id);

  /// Euler properties:
  insert<Temperature>         (p, "T0"           , 273.15 * kelvin);
  insert<Density>             (p, "rho0"         , 1.225 * kilogram
                                                        / cubic_meter);
  insert<SpecificGasConstant> (p, "Rspecific"    , 287.058 * joule
                                                   / (kilogram * kelvin));
  insert<Dimensionless>       (p, "gamma"        , 1.4);
  insert<Dimensionless>       (p, "Minfty"       , 0.1);

  /// Navier-Stokes properties:
  insert<Dimensionless>       (p, "ReInfty"      , 300);
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

  const Ind maxNoCells =  1.2 * no_solver_cells_with_gc<nd>(minRefLevel);

  /// \todo Replace with Boost.Log
  std::cerr << "fv max no cells: " << maxNoCells
            << " #leafs: " << no_leaf_nodes<nd>(minRefLevel)
            << " #faces: " << no_faces<nd>()
            << " #nodes/face: "
            << no_nodes_per_face<nd>(minRefLevel) << "\n";

  InitialDomain initialDomain = [&](const NumA<nd>) {
    return true;
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

/// \brief Flow past a moving sphere testcase
///
/// Reynolds number: Re = 80
/// Initial condition: free-stream values
/// BCs: free-stream inflow/outflow + Neumann walls + AdiabaticNoSlip
TEST(cns_fv_solver, flow_past_sphere) {
  using namespace grid::helpers::cube;

  /// Create solver
  auto cnsSolver = CNSSolver<nd> {
    cnsSolverIdx,
    cns_properties<nd>()
  };

  /// Create Heat solver:
  auto heatSolver = HeatSolver<nd> {
    heatSolverIdx, heat_properties<nd>()
  };

  /// Make moving sphere
  const Num angular_frequency = 0.05 * math::pi;
  const Num amplitude = 1.0;
  std::function<NumA<nd>()> x_center_cns = [&]() {
    NumA<nd> tmp = NumA<nd>::Constant(0.0);
    tmp(0) = -3.;
    tmp(1) = amplitude * std::sin(angular_frequency * cnsSolver.time());
    return tmp;
  };
  std::function<Num(SInd)> velocity_cns = [&](SInd d) {
    NumA<nd> v = NumA<nd>::Constant(0.);
    v(1) = amplitude * angular_frequency * std::cos(angular_frequency
                                                    * cnsSolver.time());
    return v(d);
  };
  std::function<Num(SInd)> acceleration_cns = [&](SInd d) {
    NumA<nd> a = NumA<nd>::Constant(0.);
    a(1) = - amplitude * std::pow(angular_frequency, 2)
           * std::sin(angular_frequency * cnsSolver.time());
    return a(d);
  };
  Num radius = 1.0;

  auto moving_sphere_cns =
      geometry::make_geometry<geometry::implicit::moving::Sphere<nd>>(x_center_cns, radius);

  /// Define initial condition
  auto ic = [&](const NumA<nd> x) {
    NumA<cns_nvars> pv = NumA<cns_nvars>::Zero();
    pv(CNSV::u(0))  = cnsSolver.quantities.u_infinity();
    pv(CNSV::p())   = cnsSolver.quantities.p_infinity();
    pv(CNSV::rho()) = cnsSolver.quantities.rho_infinity();
    if (x(1) > -0.25 && x(1) < 0.25) {
      Num T_wall = T_scaling * cnsSolver.quantities.T_infinity();
      pv(CNSV::rho()) = cnsSolver.quantities.gamma()
                     * cnsSolver.quantities.p_infinity() / T_wall;
    }

    return cnsSolver.cv(pv);
  };
  cnsSolver.set_initial_condition(ic);

  /// Create far-field boundary conditions
  auto nBc = cns_physics::bc::Neumann<CNSSolver<nd>>(cnsSolver);

  /// Create inflow-outflow boundary conditions
  auto inflow_with_hot_jet = [&](const CellIdx g, const SInd v) {
    NumA<cns_nvars> pv = NumA<cns_nvars>::Zero();
    pv(CNSV::u(0))  = cnsSolver.quantities.u_infinity();
    pv(CNSV::p())   = cnsSolver.quantities.p_infinity();
    pv(CNSV::rho()) = cnsSolver.quantities.rho_infinity();
    /// rho = gamma * p / T
    NumA<nd> x = cnsSolver.cells().x_center.row(g);
    if (x(1) > -0.25 && x(1) < 0.25) {
      /// rho = gamma * p / T
      Num T_wall = T_scaling * cnsSolver.quantities.T_infinity();
      pv(CNSV::rho()) = cnsSolver.quantities.gamma()
                     * cnsSolver.quantities.p_infinity() / T_wall;
    }
    return pv(v);
  };
  auto iBc = cns_physics::bc::external::Inflow<CNSSolver<nd>>
             (cnsSolver, inflow_with_hot_jet);

  auto free_stream = [&](const CellIdx, const SInd v) {
    NumA<cns_nvars> pv = NumA<cns_nvars>::Zero();
    pv(CNSV::u(0))  = cnsSolver.quantities.u_infinity();
    pv(CNSV::p())   = cnsSolver.quantities.p_infinity();
    pv(CNSV::rho()) = cnsSolver.quantities.rho_infinity();
    return pv(v);
  };
  auto oBc = cns_physics::bc::external::Outflow<CNSSolver<nd>>
             (cnsSolver, free_stream);

  /// Free-stream boundaries:
  auto bcs = std::make_tuple(iBc, oBc, nBc, nBc);
  solver::fv::append_bcs(cnsSolver, hgrid.root_cell(), bcs);

  /// Heat Solver
  std::function<NumA<nd>()> x_center_heat = [&]() {
    const Num scaling = Tr0 * cnsSolver.quantities.Pr0()
    * cnsSolver.quantities.Re0();
    const Num heat_to_cns_time = heatSolver.time() / scaling;
    NumA<nd> tmp = NumA<nd>::Constant(0.0);
    tmp(0) = -3.;
    tmp(1) = amplitude * std::sin(angular_frequency * heat_to_cns_time);
    return tmp;
  };
  std::function<Num(SInd)> velocity_heat = [&](SInd d) {
    const Num scaling = Tr0 * cnsSolver.quantities.Pr0()
    * cnsSolver.quantities.Re0();
    const Num heat_to_cns_time = heatSolver.time() / scaling;
    NumA<nd> v = NumA<nd>::Constant(0.);
    v(1) = amplitude * angular_frequency * std::cos(angular_frequency
                                                    * heat_to_cns_time);
    return v(d);
  };

  auto moving_sphere_heat =
      geometry::make_geometry<
        geometry::implicit::adaptors::Invert<
          geometry::implicit::moving::Sphere<nd>
        >
      >(x_center_heat, radius);

  /// Heat initial condition:
  auto constant_ic = [&](const NumA<nd>) {
    return NumA<heat_nvars>::Constant(cnsSolver.quantities.T_infinity());
  };
  heatSolver.set_initial_condition(constant_ic);

  /// Coupling moving boundary conditions  \todo
  auto mBh = heat_physics::bc::coupling::MBCNS<HeatSolver<nd>, CNSSolver<nd>>{
    heatSolver, cnsSolver, moving_sphere_heat, velocity_heat
  };
  heatSolver.append_bc(solver::fv::bc::Interface<nd> {
      "movingBoundary",
      moving_sphere_heat,
      heatSolver,
      mBh
  });

  /// Coupling Moving wall boundary condition: \todo
  auto mBcns = cns_physics::bc::coupling::MBHeat<CNSSolver<nd>,HeatSolver<nd>>{
    cnsSolver, heatSolver, moving_sphere_cns, velocity_cns, acceleration_cns
  };
  cnsSolver.append_bc(solver::fv::bc::Interface<nd> {
      "movingHeatCoupling",
      moving_sphere_cns,
      cnsSolver,
      mBcns
  });

  /// Initialize:
  initialize(hgrid, heatSolver, cnsSolver);
  heatSolver.impose_initial_condition_on_all_cells();
  cnsSolver.impose_initial_condition_on_all_cells();
  write_domains(hgrid, heatSolver, cnsSolver);

  /// Main loop:
  while (!solver_finished(maxNoTimeSteps, cnsSolver, heatSolver)) {
    solver::fv::coupling::apply(Tr0, cnsSolver, heatSolver);
    write_timestep(cnsSolver, heatSolver);

    cnsSolver.solve();
    heatSolver.solve();

    write_output(outputInterval, cnsSolver, heatSolver);

    if (solution_diverged(cnsSolver, heatSolver)) {
      write_domains(cnsSolver, heatSolver);
      TERMINATE("Solution diverged!");
    }
  }
  write_domains(cnsSolver, heatSolver);
}

////////////////////////////////////////////////////////////////////////////////
#undef ENABLE_DBG_
////////////////////////////////////////////////////////////////////////////////
