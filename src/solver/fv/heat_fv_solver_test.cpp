/// Includes:
#include "heat.hpp"
#include "../../geometry/geometry.hpp"
/// External Includes:
#include "gtest/gtest.h"
/// Options:
#define ENABLE_DBG_ 0
#include "../../misc/dbg.hpp"
////////////////////////////////////////////////////////////////////////////////

static const SInd nd = 3;

namespace physics = solver::fv::heat;
template<SInd nd, class S> using Physics = physics::Physics<nd,S>;
using Solver = solver::fv::Solver<nd,Physics, physics::flux::three_point>;
static const SInd nvars = Solver::nvars;
using InitialCondition = Solver::InitialCondition;
using V = typename Solver::V;
static const fvSolverId = SolverIdx{0};
const grid::RootCell<nd> rootCell(NumA<nd>::Constant(0),NumA<nd>::Constant(1));

const SInd minRefLevel = 6;
Grid<nd> test_grid(grid::helpers::cube::properties<nd>(rootCell,minRefLevel));

io::Properties fv_properties(Grid<nd>* grid, InitialCondition initialCondition,
                             const Num timeEnd) {
  using namespace grid::helpers::cube;
  const Ind maxNoCells = no_solver_cells_with_gc<nd>(minRefLevel);

  // Some log:
  std::cerr << "fv max no cells: " << maxNoCells
            << " #leafs: " << no_leaf_nodes<nd>(minRefLevel)
            << " #faces: " << no_faces<nd>()
            << " #nodes/face: "
            << no_nodes_per_face<nd>(minRefLevel) << "\n";

  io::Properties properties;
  io::insert_property<Grid<nd>*>(properties,"grid", grid);
  io::insert_property<Ind>(properties,"maxNoCells", maxNoCells);
  io::insert_property<bool>(properties,"restart",false);
  io::insert_property<InitialCondition>(properties,"initialCondition",initialCondition);
  io::insert_property<Num>(properties,"timeEnd",timeEnd);
  io::insert_property<Num>(properties,"conductivity",1.0);
  io::insert_property<Num>(properties,"CFL",0.2);
  return properties;
}

template<class BC>
io::Properties make_bcs(Solver& fvSolver, BC&& bc) {
  using namespace grid::helpers; using namespace cube; using namespace edge;
  SInd solverId = fvSolver.solver_id();
  auto boundaries
      =  make_boundaries<nd>()(solverId, rootCell, make_conditions<nd>()(std::forward<BC>(bc)));
  io::Properties properties;
  io::insert_property<decltype(boundaries)>(properties,"boundaries",boundaries);
  return properties;
}

void run_sim(Grid<nd>& grid, Solver& fvSolver, io::Properties boundaries) {

  grid.read_boundaries(boundaries);

  grid.initialize();
  write_domain("grid",&grid);

  fvSolver.initialize();
  write_domain(&fvSolver);

  while(fvSolver.time() < fvSolver.final_time() && fvSolver.step() < 1000) {
  std::cerr << "step: " << fvSolver.step() << " dt: " << fvSolver.dt()
            << " time: " << fvSolver.time() << "\n";
  fvSolver.solve();
  if(fvSolver.step() % 10 == 0 /*|| fvSolver.step() > 23*/) {
     write_domain(&fvSolver);
    }
   }
  write_domain(&fvSolver);
}

template<class BC>
void run_sim(Grid<nd>& grid, Solver& fvSolver, BC&& bc) {
  run_sim(grid,fvSolver,make_bcs(fvSolver,bc));
}

// TEST(heat_fv_solver, constant_ic_dirichlet_bc) {

//   auto constant_ic = [](const NumA<nd>) {
//     NumA<nvars> vars = NumA<nvars>::Zero();
//     return vars;
//   };

//   Solver fvSolver(fvSolverId,fv_properties(&test_grid,constant_ic,0.25));
//   auto bc = boundary::make_condition<physics::bc::Dirichlet<Solver>>(fvSolver,[](Ind,SInd){return 1.0;});
//   run_sim(test_grid,fvSolver,bc);
// }

TEST(heat_fv_solver, one_dimensional_dirichlet) {
  using DBc = physics::bc::Dirichlet<Solver>;
  using NBc = physics::bc::Neumann<Solver>;
  auto constant_ic = [](const NumA<nd>) { return NumA<nvars>::Zero(); };

  Solver fvSolver(fvSolverId,fv_properties(&test_grid,constant_ic,0.25));

  auto dBc1 = boundary::make_condition<DBc>(fvSolver,[](Ind,SInd){return 1.0;});
  auto dBc2 = boundary::make_condition<DBc>(fvSolver,[](Ind,SInd){return 2.0;});
  auto dBc3 = boundary::make_condition<DBc>(fvSolver,[](Ind,SInd){return 3.0;});
  auto dBc4 = boundary::make_condition<DBc>(fvSolver,[](Ind,SInd){return 4.0;});
  auto dBc5 = boundary::make_condition<DBc>(fvSolver,[](Ind,SInd){return 5.0;});
  auto dBc6 = boundary::make_condition<DBc>(fvSolver,[](Ind,SInd){return 6.0;});;
  auto nBc = boundary::make_condition<NBc>(fvSolver,[](Ind,SInd){return 0.0;});

  auto bcs = std::make_tuple(dBc1,dBc1,nBc,nBc,nBc,nBc);
  //auto bcs = std::make_tuple(dBc,dBc,dBc,dBc,dBc,dBc);
  auto bcs = std::make_tuple(dBc1,dBc2,dBc3,dBc4,dBc5,dBc6);

  SInd solverId = fvSolver.solver_id();
  auto boundaries = grid::helpers::cube::make_boundaries<nd>()(solverId, rootCell, bcs);
  io::Properties properties;
  io::insert_property<decltype(boundaries)>(properties,"boundaries",boundaries);


  run_sim(test_grid,fvSolver,properties);
}

////////////////////////////////////////////////////////////////////////////////
#undef ENABLE_DBG_
////////////////////////////////////////////////////////////////////////////////
