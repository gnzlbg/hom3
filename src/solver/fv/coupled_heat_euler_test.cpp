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

static const SInd nd = 3;
static const auto eulerSolverId = SolverIdx{0};
static const auto heatSolverId = SolverIdx{1};
static const Num timeEnd = 0.25;
static const Num cfl = 0.2;

// Grid:
const grid::RootCell<nd> rootCell(NumA<nd>::Constant(0),NumA<nd>::Constant(1));
const SInd minRefLevel = 4;

auto cubeH = geometry::make_geometry<geometry::implicit::Square<nd>>(NumA<nd>{0.5,0.5,0.5},
                                                                    NumA<nd>{0.2,0.2,0.2});
auto cubeE = geometry::make_geometry<geometry::implicit::Square<nd>>(NumA<nd>{0.5,0.5,0.5},
                                                                     NumA<nd>{0.26,0.26,0.26});
Grid<nd> test_grid(grid::helpers::cube::properties<nd>(rootCell,minRefLevel,2));

// Euler Fv Solver:
namespace euler_physics = solver::fv::euler;
template<SInd nd, class S> using EulerPhysics = euler_physics::Physics<nd,S>;
using EulerSolver = solver::fv::Solver<nd,EulerPhysics, euler_physics::flux::ausm>;

io::Properties euler_properties(Grid<nd>* grid) {
  using namespace grid::helpers::cube;
  const Ind maxNoCells = no_solver_cells_with_gc<nd>(minRefLevel)*2;
  static const SInd nvars = EulerSolver::nvars;
  using V = euler_physics::Indices<nd>;

  EulerSolver::InitialDomain initialDomain = [=](const NumA<nd> x){
    return (*cubeH)(x) > 0 ? true : false;
  };

  EulerSolver::InitialCondition initialCondition = [](const NumA<nd>) {
    NumA<nvars> pvars = NumA<nvars>::Zero();
    pvars(V::rho()) = 1.0;
    for(SInd d = 0; d < nd; ++d) {
      pvars(V::u(d)) = 0.0;
    }
    pvars(V::p()) = 1.0;
    return euler_physics::cv<nd>(pvars,1.4);
  };

  io::Properties properties;
  io::insert_property<Grid<nd>*>(properties,"grid", grid);
  io::insert_property<Ind>(properties,"maxNoCells", maxNoCells);
  io::insert_property<bool>(properties,"restart",false);
  io::insert_property<EulerSolver::InitialDomain>(properties,"initialDomain",initialDomain);
  io::insert_property<EulerSolver::InitialCondition>(properties,"initialCondition",initialCondition);
  io::insert_property<Num>(properties,"timeEnd",timeEnd);
  io::insert_property<Num>(properties,"gamma",1.4);
  io::insert_property<Num>(properties,"CFL",cfl);
  return properties;
}

Grid<nd>::Boundaries euler_bcs(EulerSolver& solver) {
  using namespace grid::helpers; using namespace cube; using namespace edge;
  auto bc
  = boundary::make_condition<euler_physics::bc::Neumann<EulerSolver>>(solver);
  auto bcs = make_conditions<nd>()(bc);
  auto boundaries =  make_boundaries<nd>()(eulerSolverId, rootCell, bcs);
  return boundaries;
}

// Heat Fv Solver:
namespace heat_physics = solver::fv::heat;
template<SInd nd, class S> using HeatPhysics = heat_physics::Physics<nd,S>;
using HeatSolver = solver::fv::Solver<nd,HeatPhysics, heat_physics::flux::three_point>;

io::Properties heat_properties(Grid<nd>* grid) {
  using namespace grid::helpers::cube;
  const Ind maxNoCells = no_solver_cells_with_gc<nd>(minRefLevel);
  static const SInd nvars = HeatSolver::nvars;
  using V = heat_physics::Indices<nd>;

  HeatSolver::InitialDomain initialDomain = [=](const NumA<nd> x){
    return math::apprx((*cubeH)(x),0.) || (*cubeH)(x) < 0. ? true : false;
  };

  HeatSolver::InitialCondition initialCondition
      = [](const NumA<nd>) { return NumA<nvars>::Constant(1.0); };

  io::Properties properties;
  io::insert_property<Grid<nd>*>(properties,"grid", grid);
  io::insert_property<Ind>(properties,"maxNoCells", maxNoCells);
  io::insert_property<bool>(properties,"restart",false);
  io::insert_property<HeatSolver::InitialDomain>(properties,"initialDomain",initialDomain);
  io::insert_property<HeatSolver::InitialCondition>(properties,"initialCondition",initialCondition);
  io::insert_property<Num>(properties,"timeEnd",timeEnd);
  io::insert_property<Num>(properties,"conductivity",1.0);
  io::insert_property<Num>(properties,"CFL",cfl);
  return properties;
}

int main() {

  EulerSolver eulerSolver(eulerSolverId,euler_properties(&test_grid));
  HeatSolver heatSolver(heatSolverId,heat_properties(&test_grid));

  auto boundaries = euler_bcs(eulerSolver);

  // Euler-Heat interface conditions
  //auto ehBc =  boundary::make_condition<solver::fv::euler::bc::Neumann<EulerSolver>>
  //             (eulerSolver, [](Ind,SInd)->Num{ return 0; });
   auto ehBc = boundary::make_condition<solver::fv::coupling::bc::EulerHeatWall<EulerSolver,HeatSolver>>
               (eulerSolver,heatSolver);
  // auto heBc = boundary::make_condition<solver::fv::heat::bc::Dirichlet<HeatSolver>>
  //             (heatSolver,[](Ind,SInd){return 1.0;});

  auto heBc = boundary::make_condition<solver::fv::coupling::bc::HeatEulerWall<HeatSolver,EulerSolver>>
               (heatSolver,eulerSolver);

  auto eulerCouplingBc = boundary::Interface<nd>{"eulerCoupling",cubeE,{eulerSolverId},{ehBc}};
  auto heatCouplingBc = boundary::Interface<nd>{"heatCoupling",geometry::implicit::adaptors::invert2(cubeH),{heatSolverId},{heBc}};
  //auto boundaries = std::vector<boundary::Interface<nd>>{ heatCouplingBc };
  boundaries.push_back(eulerCouplingBc);
  boundaries.push_back(heatCouplingBc);


  io::Properties boundaryProperties;
  io::insert_property<decltype(boundaries)>(boundaryProperties,"boundaries",boundaries);

  test_grid.read_boundaries(boundaryProperties);

  test_grid.initialize();
  write_domain("grid",&test_grid);

  eulerSolver.initialize();
  heatSolver.initialize();

  write_domain(&eulerSolver);
  write_domain(&heatSolver);

  auto print_total_energy = [&](){
    Num eulerE0 = 0;
    for(auto&& cId : eulerSolver.internal_cells()) {
      eulerE0 += eulerSolver.cells().lhs(cId,EulerSolver::V::E());
    }
    Num heatE0 = 0;
    for(auto&& cId : heatSolver.internal_cells()) {
      heatE0 += heatSolver.cells().lhs(cId,HeatSolver::V::T());
    }
    std::cerr << "E_euler = " << eulerE0
    << " | E_heat = " << heatE0
    << " | E_total = " << eulerE0 + heatE0
    << std::endl;
  };

  print_total_energy();
  
  while(
        eulerSolver.time() < eulerSolver.final_time()
          &&
        heatSolver.time() < heatSolver.final_time()
        ) {

    //const Num dt = eulerSolver.min_dt();
    //const Num dt = heatSolver.min_dt();
    const Num dt = std::min(eulerSolver.min_dt(),heatSolver.min_dt());
    eulerSolver.dt(dt);
    heatSolver.dt(dt);

    std::cerr
        << "t_euler = " << eulerSolver.time()
        << " | t_heat = " << heatSolver.time()
        << " | dt = " << dt;

    eulerSolver.solve();
    heatSolver.solve();

    print_total_energy();

    if(//eulerSolver.step() % 1 == 0
           heatSolver.step() % 1 == 0
       ) {
      write_domain(&eulerSolver);
      write_domain(&heatSolver);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
#undef ENABLE_DBG_
////////////////////////////////////////////////////////////////////////////////
