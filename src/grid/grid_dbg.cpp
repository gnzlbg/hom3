////////////////////////////////////////////////////////////////////////////////
/// Library Includes:
#include <algorithm>
#include <vector>
/// External Includes:
#include "globals.hpp"
//#include "../cell/testgridcell.hpp"
#include "containers/hierarchical.hpp"
#include "grid.hpp"
#include "geometry/implicit.hpp"

// EXPECT_DEATH(cells[-1].myNum() = -2.0, "[\\S\\s]+");

template<class Grid>
void write_grid(const std::string fName, Grid* grid) {

  static const int nd = Grid::no_dimensions();

  io::Vtk<3,io::format::binary> out(grid, fName, io::precision::standard());

  out << io::stream("cellIds",1,[](const Ind cId, const SInd){ return cId; });

  NumA<nd> center = { 0.5, 0.5, 0.5 };
  auto sphere = geometry::implicit::Sphere<nd> (center, 0.25);
  out << io::stream("levelSet",1,[&](const Ind cId, const SInd){
      return sphere(grid->cell_coordinates(cId));
    });
  out << io::stream("nghbrs",grid->nodes().no_samelvl_nghbr_pos(),
                    [&](const Ind cId, const SInd pos){
                      return grid->nodes().find_samelvl_nghbr(cId,pos);
                    });
  out << io::stream("localIds",grid->nodes().container_capacity(),
                    [&](const Ind cId, const SInd c){
                      return grid->nodes().local_id(cId,c);
                    });

}

int main() {

  //typedef Grid<2,1,grid::generation::MinLevel> TestGrid;
  typedef Grid<3> TestGrid;
  typedef grid::generation::MinLevel MeshGen;

  //Ind maxNoGridNodes = 600;
  constexpr Ind minRefinementLevel = 7;
  constexpr Ind maxNoGridNodes = grid::helpers::unit_cube::no_nodes<3>(minRefinementLevel);
  constexpr SInd maxNoContainers = 1;
  // NumA<2> min = {0,0};
  // NumA<2> max = {1,1};
  // grid::RootCell<2> rootCell(min,max);
  const NumA<3> min = {0,0,0};
  const NumA<3> max = {1,1,1};
  const grid::RootCell<3> rootCell(min,max);

  io::Properties meshGeneration;
  io::insert_property<Ind>(meshGeneration,"level",minRefinementLevel);

  MeshGen mesh_gen(meshGeneration);

  io::Properties properties;
  io::insert_property<grid::RootCell<3>>(properties,"rootCell",rootCell);
  io::insert_property<Ind>(properties,"maxNoGridNodes",maxNoGridNodes);
  io::insert_property<std::function<void(TestGrid*)>>(properties,"meshGeneration",mesh_gen);
  io::insert_property<SInd>(properties,"maxNoContainers",maxNoContainers);

  TestGrid grid(properties, grid::initialize);
  std::cout << grid.nodes().size() << std::endl;
  write_grid("test_grid",&grid);

}
