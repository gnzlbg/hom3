/// Includes:
#include "heat.hpp"
#include "euler.hpp"
#include "coupling.hpp"
#include "../../grid/helpers.hpp"
#include "../../geometry/geometry.hpp"
/// External Includes:
#include "gtest/gtest.h"
/// Options:
#define ENABLE_DBG_ 0
#include "../../misc/dbg.hpp"
////////////////////////////////////////////////////////////////////////////////

using namespace hom3;

static const SInd nd = 3;
static const auto euler0SolverId = SolverIdx{0};
static const auto euler1SolverId = SolverIdx{1};
static const Num timeEnd = 0.25;
static const Num cfl = 0.2;

// Grid:
const grid::RootCell<nd> rootCell(NumA<nd>::Constant(0),NumA<nd>::Constant(1));
const SInd minRefLevel = 4;

auto cubeH = geometry::make_geometry<geometry::implicit::Square<nd>>(NumA<nd>{0.5,0.5,0.5},
                                                                    NumA<nd>{0.2,0.2,0.2});
auto cubeE = geometry::make_geometry<geometry::implicit::Square<nd>>(NumA<nd>{0.5,0.5,0.5},
                                                                     NumA<nd>{0.26,0.26,0.26});
grid::Grid<nd> test_grid(grid::helpers::cube::properties<nd>(rootCell,minRefLevel,2));

// EulerG Fv Solver:
namespace physics = solver::fv::euler;
template<SInd nd, class S> using EulerPhysics = physics::Physics<nd,S>;
using EulerSolver = solver::fv::Solver<nd,EulerPhysics, physics::flux::ausm>;

io::Properties euler_properties(grid::Grid<nd>* grid, SolverIdx solverId) {
  using namespace grid::helpers::cube;
  const Ind maxNoCells = no_solver_cells_with_gc<nd>(minRefLevel)*2;
  static const SInd nvars = EulerSolver::nvars;
  using V = physics::Indices<nd>;

  EulerSolver::InitialDomain initialDomain;

  if(solverId() == 0) {
    initialDomain = [=](const NumA<nd> x){
      return (*cubeH)(x) > 0 ? true : false;
    };
  } else {
    initialDomain = [=](const NumA<nd> x){
      return math::apprx((*cubeH)(x),0.) || (*cubeH)(x) < 0. ? true : false;
    };
  }

  static constexpr SInd dir = 0;
  static constexpr Num angle = 0;

  EulerSolver::InitialCondition initialCondition = [=](const NumA<nd> x) {
    const Num alpha = angle * (2 * math::pi) / 360;
    NumA<nd> x_rot;
    x_rot(0) = std::cos(alpha) * x(0) - std::sin(alpha) * x(1);
    x_rot(1) = std::sin(alpha) * x(0) + std::cos(alpha) * x(1);
    NumA<nvars> pvars;

    const Num pos = math::apprx<Num>(alpha,0) ? 0.5 : std::sqrt(2)/2;

    if(x_rot(dir) > pos) {
      pvars(V::rho()) = 1.0;
      for(SInd d = 0; d < nd; ++d) {
        pvars(V::u(d)) = 0;
      }
      pvars(V::p()) = 1.0;
    } else {
      pvars(V::rho()) = 0.125;
      for(SInd d = 0; d < nd; ++d) {
        pvars(V::u(d)) = 0;
      }
      pvars(V::p()) = 0.1;
    }
    return physics::cv<nd>(pvars,1.4);
  };

  io::Properties properties;
  io::insert_property<grid::Grid<nd>*>(properties,"grid", grid);
  io::insert_property<Ind>(properties,"maxNoCells", maxNoCells);
  io::insert_property<bool>(properties,"restart",false);
  io::insert_property<EulerSolver::InitialDomain>(properties,"initialDomain",initialDomain);
  io::insert_property<EulerSolver::InitialCondition>(properties,"initialCondition",initialCondition);
  io::insert_property<Num>(properties,"timeEnd",timeEnd);
  io::insert_property<Num>(properties,"gamma",1.4);
  io::insert_property<Num>(properties,"CFL",cfl);
  return properties;
}



grid::Grid<nd>::Boundaries euler0_bcs(EulerSolver& solver) {
  using namespace grid::helpers; using namespace cube; using namespace edge;
  auto bc
  = boundary::make_condition<physics::bc::Neumann<EulerSolver>>(solver);
  auto bcs = make_conditions<nd>()(bc);
  auto boundaries =  make_boundaries<nd>()(euler0SolverId, rootCell, bcs);
  return boundaries;
}

int main() {

  EulerSolver euler0Solver(euler0SolverId,euler_properties(&test_grid,euler0SolverId));
  EulerSolver euler1Solver(euler1SolverId,euler_properties(&test_grid,euler1SolverId));


  auto boundaries = euler0_bcs(euler0Solver);

  auto e0e1Bc = boundary::make_condition<solver::fv::coupling::bc::EulerEulerWall<EulerSolver>>
                (euler0Solver,euler1Solver);
  auto e1e0Bc = boundary::make_condition<solver::fv::coupling::bc::EulerEulerWall<EulerSolver>>
                (euler1Solver,euler0Solver);

  auto euler0CouplingBc = boundary::Interface<nd>{"euler0Coupling",cubeE,{euler0SolverId},{e0e1Bc}};
  auto euler1CouplingBc = boundary::Interface<nd>{"euler1Coupling",geometry::implicit::adaptors::invert2(cubeH),{euler1SolverId},{e1e0Bc}};

  boundaries.push_back(euler0CouplingBc);
  boundaries.push_back(euler1CouplingBc);

  io::Properties boundaryProperties;
  io::insert_property<decltype(boundaries)>(boundaryProperties,"boundaries",boundaries);

  test_grid.read_boundaries(boundaryProperties);

  test_grid.initialize();

  euler0Solver.initialize();
  euler1Solver.initialize();

  write_domain("grid",&test_grid);

  write_domain(&euler0Solver);
  write_domain(&euler1Solver);

  while(
        euler0Solver.time() < euler0Solver.final_time()
        &&
        euler1Solver.time() < euler1Solver.final_time()
        ) {

    const Num dt = std::min(euler0Solver.min_dt(),euler1Solver.min_dt());
    euler0Solver.dt(dt);
    euler1Solver.dt(dt);

    std::cerr
        << "t_euler0 = " << euler0Solver.time()
        << " | t_euler1 = " << euler1Solver.time()
        << " | dt = " << dt << std::endl;

    euler0Solver.solve();
    euler1Solver.solve();

    if(euler0Solver.step() % 1 == 0
       ||
       euler1Solver.step() % 1 == 0
       ) {
      write_domain(&euler0Solver);
      write_domain(&euler1Solver);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
#undef ENABLE_DBG_
////////////////////////////////////////////////////////////////////////////////
