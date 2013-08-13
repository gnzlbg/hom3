/// Includes:
#include "euler.hpp"
#include "../../geometry/geometry.hpp"
/// External Includes:
#include "gtest/gtest.h"
/// Options:
#define ENABLE_DBG_ 0
#include "../../misc/dbg.hpp"
////////////////////////////////////////////////////////////////////////////////

using namespace hom3;

static const SInd nd = 3;

namespace physics = solver::fv::euler;
template<SInd nd, class S> using Physics = physics::Physics<nd,S>;
using Solver = solver::fv::Solver<nd, Physics, physics::flux::ausm>;
static const SInd nvars = Solver::nvars;
using InitialCondition = Solver::InitialCondition;
using V = typename Solver::V;
static const SInd fvSolverId = 0;
const grid::RootCell<nd> rootCell(NumA<nd>::Constant(0),NumA<nd>::Constant(1));

const SInd minRefLevel = 6;
Grid<nd> test_grid(grid::helpers::cube::properties<nd>(rootCell,minRefLevel));

template<SInd nd, class InitCon>
io::Properties fv_properties(Grid<nd>* grid, InitCon initCon,
                             const Num timeEnd) {
  using namespace grid::helpers::cube;
  using InitialCondition = typename solver::fv::Solver<nd,Physics, physics::flux::ausm>::InitialCondition;
  const Ind maxNoCells = no_solver_cells_with_gc<nd>(minRefLevel);
  using InitialDomain = std::function<bool(const NumA<nd>)>;
  // Some log:
  std::cerr << "fv max no cells: " << maxNoCells
            << " #leafs: " << no_leaf_nodes<nd>(minRefLevel)
            << " #faces: " << no_faces<nd>()
            << " #nodes/face: "
            << no_nodes_per_face<nd>(minRefLevel) << "\n";

  InitialCondition initialCondition = initCon;
  InitialDomain initialDomain = [](const NumA<nd>)->bool{return true;};

  io::Properties properties;
  io::insert_property<Grid<nd>*>(properties,"grid", grid);
  io::insert_property<Ind>(properties,"maxNoCells", maxNoCells);
  io::insert_property<bool>(properties,"restart",false);
  io::insert_property<InitialDomain>(properties,"initialDomain",initialDomain);
  io::insert_property<InitialCondition>(properties,"initialCondition",initialCondition);
  io::insert_property<Num>(properties,"timeEnd",timeEnd);
  io::insert_property<Num>(properties,"gamma",1.4);
  io::insert_property<Num>(properties,"CFL",0.2);
  return properties;
}

template<class Solver, class BC>
io::Properties make_bcs(Solver& fvSolver, BC&& bc) {
  using namespace grid::helpers; using namespace cube; using namespace edge;
  SInd solverId = fvSolver.solver_id();
  auto boundaries
      =  make_boundaries<nd>()(solverId, rootCell, make_conditions<nd>()(std::forward<BC>(bc)));
  io::Properties properties;
  io::insert_property<decltype(boundaries)>(properties,"boundaries",boundaries);
  return properties;
}

template<SInd nd, class Solver, class BC>
void run_sim(Grid<nd>& grid, Solver& fvSolver, BC&& bc) {

  grid.read_boundaries(make_bcs(fvSolver,bc));

  grid.initialize();
  write_domain("grid",&grid);

  fvSolver.initialize();
  write_domain(&fvSolver);

  while(fvSolver.time() < fvSolver.final_time() && fvSolver.step() < 1000) {
  std::cerr << "step: " << fvSolver.step() << " dt: " << fvSolver.dt()
            << " time: " << fvSolver.time() << "\n";
  fvSolver.solve();
  if(fvSolver.step() % 1 == 0 /*|| fvSolver.step() > 23*/) {
     write_domain(&fvSolver);
    }
   }
  write_domain(&fvSolver);
}

TEST(euler_fv_solver, constant_ic) {

  auto constant_ic = [](const NumA<nd>) {
    NumA<nvars> pvars = NumA<nvars>::Zero();
    pvars(V::rho()) = 1.0;
    for(SInd d = 0; d < nd; ++d) {
      pvars(V::u(d)) = 0.5;
    }
    pvars(V::p()) = 1.0;
    return physics::cv<nd>(pvars,1.4);
  };

  Solver fvSolver(fvSolverId,fv_properties<nd>(&test_grid,constant_ic,0.25));
  auto bc = boundary::make_condition<physics::bc::Neumann<Solver>>(fvSolver);
  run_sim(test_grid,fvSolver,bc);
}


TEST(euler_fv_solver, sod_shock_tube_ic) {

  const SInd dir = 0;
  const Num angle = 0;

  auto sod_shock_tube_ic = [=](const NumA<nd> x) {
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

  Solver fvSolver(fvSolverId,fv_properties<nd>(&test_grid,sod_shock_tube_ic,0.25));
  auto bc = boundary::make_condition<physics::bc::Neumann<Solver>>(fvSolver);
  run_sim(test_grid,fvSolver,bc);
}

TEST(euler_fv_solver, explosion_ic) {
  const Num radius = 0.2;
  auto explosion_ic = [=](const NumA<nd> x) {
    const Num R = radius;
    NumA<nd> xE;
    xE(0) = 0.5; xE(1) = 0.5;
      if(nd == 3) { xE(2) = 0.5; }

    NumA<nvars> pvars;
    if((x - xE).norm() < R) { // inside
      pvars(V::rho()) = 1.0;
      for(SInd d = 0; d < nd; ++d) {
        pvars(V::u(d)) = 0;
      }
      pvars(V::p()) = 1.0;
      } else { // outside
      pvars(V::rho()) = 0.125;
      for(SInd d = 0; d < nd; ++d) {
        pvars(V::u(d)) = 0;
      }
      pvars(V::p()) = 0.1;
    }
    return physics::cv<nd>(pvars,1.4);
  };

  Solver fvSolver(fvSolverId,fv_properties<nd>(&test_grid,explosion_ic,0.25));
  auto bc = boundary::make_condition<physics::bc::Neumann<Solver>>(fvSolver);
  run_sim(test_grid,fvSolver,bc);
}

TEST(euler_fv_solver, isentropic_vortex_ic) {

  using Solver2D = solver::fv::Solver<2,Physics, physics::flux::ausm>;
  const grid::RootCell<2> rootCell_2d(NumA<2>::Constant(0), NumA<2>::Constant(1));
  Grid<2> test_grid_2d(grid::helpers::cube::properties<2>(rootCell_2d,minRefLevel));

  auto isentropic_vortex_ic = [](const NumA<2> x) {
    return physics::isentropic_vortex<2>(x,0);
  };

  Solver2D fvSolver(fvSolverId,fv_properties<2>(&test_grid_2d,isentropic_vortex_ic,0.25));
  auto bc = boundary::make_condition<physics::bc::IsentropicVortex<Solver2D>>(fvSolver);
  run_sim(test_grid_2d,fvSolver,bc);
}


////////////////////////////////////////////////////////////////////////////////
#undef ENABLE_DBG_
////////////////////////////////////////////////////////////////////////////////
