////////////////////////////////////////////////////////////////////////////////
/// Library Includes:
#include <algorithm>
#include <vector>
/// External Includes:
#include "../globals.hpp"
//#include "../cell/testgridcell.hpp"
#include "../cell/container.hpp"
#include "../grid/grid.hpp"
#include "implicit_functions.hpp"

template<class Grid>
void write_grid_with_geom(const std::string fName,Grid* grid) {
  static const SInd nd = Grid::no_dimensions();

  io::Vtk out(grid, fName + "_" + std::to_string(nd), io::format::ascii(), io::precision::standard());

  out << io::stream("cellIds",1,[](const Ind cId, const SInd){ return cId; });

  const NumA<nd> center = nd == 2? { 0.5, 0.5 } : { 0.5, 0.5, 0.5 };
  const Num radius = 0.25;
  out << io::stream("levelSet",1,[&](const Ind cId, const SInd){
      const NumA<nd> x_cell = grid->cell_coordinates(cId);
      return level_set_circle(x_cell,center,radius);
    });
  out << io::stream("nghbrs",grid->nodes_.no_samelvl_nghbr_pos(),[&](const Ind cId, const SInd pos){
      return grid->nodes_.samelvl_nghbr(cId,pos);
    });
  out << io::stream("localIds",grid->nodes_.container_capacity(),[&](const Ind cId, const SInd c){
      return grid->nodes_.local_id(cId,c);
    });

  auto s = geometry::make_geometry<geometry::implicit::Sphere<nd>>(center, radius);
  out << io::stream("sphere",1,[&](const Ind cId, const SInd){
      return (*s)(grid->cell_coordinates(cId));
    });
  NumA<nd> normal = nd == 2 ? {1.0, 0.0} : {1.0, 0.0, 0.0};
  auto surface = geometry::make_geometry<
      geometry::implicit::adaptors::Invert<
                        geometry::implicit::Edge<nd>
                        >
      >(center, normal);
  out << io::stream("surface",1,[&](const Ind cId, const SInd){
      return (*surface)(grid->cell_coordinates(cId));
    });

  NumA<nd> lengths = nd == 2 ? {0.2, 0.2} : {0.2, 0.2, 0.2};
  NumA<nd> center1 = nd == 2 ? {0.25, 0.25} : {0.25, 0.25, 0.25};

  auto sq1 = geometry::make_geometry<geometry::implicit::Square<nd>>(center1, lengths);
  out << io::stream("square",1,[&](const Ind cId, const SInd){
      return (*sq1)(grid->cell_coordinates(cId));
    });

  NumA<nd> center2 = nd == 2 ? {0.75, 0.75} : {0.75, 0.75, 0.75};
  auto sq2 = geometry::make_geometry<geometry::implicit::Square<nd>>(center2, lengths);
  auto un = geometry::make_geometry
            <decltype(geometry::implicit::adaptors::make_union(sq1,sq2))>
            (geometry::implicit::adaptors::make_union(sq1,sq2));
  out << io::stream("union",1,[&](const Ind cId, const SInd){
      return (*un)(grid->cell_coordinates(cId));
    });

  for(const auto& i : grid->boundaries) {
    out << io::stream(i.name(),1,[&](const Ind cId, const SInd) {
        return i.signed_distance(grid->cell_coordinates(cId));
      });
  }
}


void test_3d() {
  static const SInd nd = 3;
  typedef Grid<3> TestGrid;
  typedef grid::generation::MinLevel MeshGen;

  Ind maxNoGridNodes = 50000;
  SInd maxNoContainers = 1;
  const NumA<3> min = {0,0,0};
  const NumA<3> max = {1,1,1};
  const grid::RootCell<3> rootCell(min,max);
  std::vector<TestGrid::Boundary> boundaries;

  const NumA<nd> center = nd == 2 ? { 0.5, 0.5 } : { 0.5, 0.5, 0.5 };
  const NumA<nd> normal = nd == 2 ? {1.0, 0.0} : {1.0, 0.0, 0.0};
  auto surface = geometry::make_geometry<geometry::implicit::Edge<nd>>(center, normal);
  boundaries.push_back(TestGrid::BoundaryInterface("bc0",surface,1));
  Ind minRefinementLevel = 4;

  io::Properties meshGeneration;
  io::insert_property<Ind>(meshGeneration,"level",minRefinementLevel);

  MeshGen mesh_gen(meshGeneration);

  io::Properties properties;
  io::insert_properties<grid::RootCell<3>>(properties,"rootCell", rootCell);
  io::insert_properties<Ind>(properties,"maxNoGridNodes",maxNoGridNodes);
  io::insert_properties<std::function<void(TestGrid*)>>(properties,"meshGeneration",mesh_gen);
  io::insert_properties<SInd>(properties,"maxNoContainers",maxNoContainers);
  io::insert_properties<std::vector<Boundary>(properties,"boundaries",boundaries);

  TestGrid grid(properties);

  write_grid_with_geom("test_grid",&grid);
}


void test_2d() {
  static const SInd nd = 2;
  typedef Grid<2> TestGrid;
  typedef grid::generation::MinLevel MeshGen;

  Ind maxNoGridNodes = 50000;
  SInd maxNoContainers = 1;
  const NumA<2> min = {0,0};
  const NumA<2> max = {1,1};
  const grid::RootCell<2> rootCell(min,max);
  Ind minRefinementLevel = 6;

  io::properties meshGeneration;
  meshGeneration.insert(io::property("level",minRefinementLevel));

  MeshGen mesh_gen(meshGeneration);
  std::vector<TestGrid::BoundaryInterface> boundaries;

   NumA<nd> center;
  switch(nd) {
    case 2: { center << 0.5, 0.5; break; }
    case 3: { center << 0.5, 0.5, 0.5; break; }
  }
  NumA<nd> normal;
  switch(nd) {
    case 2: { normal << 1.0, 0.0; break; }
    case 3: { normal << 1.0, 0.0, 0.0; break; }
    default: { TERMINATE("Wrong #of dimensions"); }
  }
  auto surface = geometry::make_geometry<geometry::implicit::Edge<nd>>(center, normal);
  boundaries.push_back(TestGrid::BoundaryInterface("bc0",surface,1));

  io::properties properties;
  properties.insert(io::property("rootCell"      , rootCell           ));
  properties.insert(io::property("maxNoGridNodes", Ind(maxNoGridNodes)));
  properties.insert(io::property("meshGeneration", std::function<void(TestGrid*)>(mesh_gen)     ));
  properties.insert(io::property("maxNoContainers", SInd(maxNoContainers)));
    properties.insert(io::property("boundaries",boundaries));

  TestGrid grid(properties);

  write_grid_with_geom("test_grid",&grid);

}

int main() {

  //typedef Grid<2,1,grid::generation::MinLevel> TestGrid;
  test_3d();
  test_2d();
}
