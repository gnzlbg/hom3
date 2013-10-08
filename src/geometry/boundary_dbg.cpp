////////////////////////////////////////////////////////////////////////////////
/// Library Includes:
#include <algorithm>
#include <vector>
/// External Includes:
#include "globals.hpp"
//#include "../cell/testgridcell.hpp"
#include "containers/sequential.hpp"
#include "grid/grid.hpp"
#include "geometry/implicit.hpp"

template<class Grid>
void write_grid_with_boundary(const std::string fName,Grid* grid) {
  static const SInd nd = Grid::no_dimensions();

  io::Vtk out(grid, fName + "_" + std::to_string(nd), io::format::ascii(), io::precision::standard());

  out << io::stream("cellIds",1,[](const Ind cId, const SInd){ return cId; });

  out << io::stream("levelSet",1,[&](const Ind cId, const SInd){
      return grid->level_set(cId);
    });
  out << io::stream("nghbrs",grid->nodes_.no_samelvl_nghbr_pos(),[&](const Ind cId, const SInd pos){
      return grid->nodes_.samelvl_nghbr(cId,pos);
    });
  out << io::stream("localIds",grid->nodes_.container_capacity(),[&](const Ind cId, const SInd c){
      return grid->nodes_.local_id(cId,c);
    });

  for(const auto& i : grid->boundaries) {
    out << io::stream(i.name() + "_ls",1,[&](const Ind cId, const SInd) {
        return i.signed_distance(grid->cell_coordinates(cId));
      });
    out << io::stream(i.name() + "_ids",1,[&](const Ind cId, const SInd) {
        return grid->is_cut_by(cId,[&](const NumA<nd> x){ return i.signed_distance(x); }) ? Ind(1) : Ind(0);
      });
  }
}

void test_2d() {
  static const SInd nd = 2;
  typedef Grid<2> TestGrid;
  typedef grid::generation::MinLevel MeshGen;

  Ind maxNoGridNodes = 50000;
  SInd maxNoContainers = 1;
  NumA<2> min = {0,0};
  NumA<2> max = {1,1};
  grid::RootCell<2> rootCell(min,max);
  Ind minRefinementLevel = 6;

  io::properties meshGeneration;
  meshGeneration.insert(io::property("level",minRefinementLevel));

  MeshGen mesh_gen(meshGeneration);
  std::vector<TestGrid::BoundaryInterface> boundaries;

  { // Bottom
    NumA<nd> center; center << 0.5, 0.00001;
    NumA<nd> normal; normal << 0.0, 1.0;
    auto surface = geometry::make_geometry<geometry::implicit::Edge<nd>> (center, normal);
    boundaries.push_back(TestGrid::BoundaryInterface("bc_bottom",surface,1));
  }

  { // Top
    NumA<nd> center; center << 0.5, 0.99999;
    NumA<nd> normal; normal << 0.0, -1.0;
    auto surface = geometry::make_geometry<geometry::implicit::Edge<nd>> (center, normal);
    boundaries.push_back(TestGrid::BoundaryInterface("bc_top",surface,1));
  }

  { // Left
    NumA<nd> center; center << 0.0001, 0.5;
    NumA<nd> normal; normal << 1.0, 0.0;
    auto surface = geometry::make_geometry<geometry::implicit::Edge<nd>> (center, normal);
    boundaries.push_back(TestGrid::BoundaryInterface("bc_left",surface,1));
  }
  
  { // Right
    NumA<nd> center; center << 0.99999, 0.5;
    NumA<nd> normal; normal << -1.0, 0.0;
    auto surface = geometry::make_geometry<geometry::implicit::Edge<nd>> (center, normal);
    boundaries.push_back(TestGrid::BoundaryInterface("bc_right",surface,1));
  }

  io::properties properties;
  properties.insert(io::property("rootCell"      , rootCell           ));
  properties.insert(io::property("maxNoGridNodes", Ind(maxNoGridNodes)));
  properties.insert(io::property("meshGeneration", std::function<void(TestGrid*)>(mesh_gen)     ));
  properties.insert(io::property("maxNoContainers", SInd(maxNoContainers)));
  properties.insert(io::property("boundaries",boundaries));

  TestGrid grid(properties);

  write_grid_with_boundary("boundary_test",&grid);

}

int main() {

  test_2d();
}
