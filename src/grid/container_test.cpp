/// Includes:
#include "grid/grid.hpp"
/// External Includes:
#include "gtest/gtest.h"
/// Options:
#define ENABLE_DBG_ 0
#include "misc/dbg.hpp"
////////////////////////////////////////////////////////////////////////////////

typedef Grid<2> Test2DGrid;
typedef Grid<3> Test3DGrid;
typedef grid::generation::MinLevel MeshGen;

io::Properties small_3D_grid() {
  Ind maxNoGridNodes = 100;
  SInd maxNoContainers = 1;
  NumA<3> min = {0,0,0};
  NumA<3> max = {1,1,1};
  grid::RootCell<3> rootCell(min,max);
  Ind minRefinementLevel = 2;
  io::Properties meshGeneration;
  io::insert_property<Ind>(meshGeneration,"level",minRefinementLevel);
  MeshGen mesh_gen(meshGeneration);
  io::Properties properties;
  io::insert_property<grid::RootCell<3>>(properties, "rootCell", rootCell);
  io::insert_property<Ind>(properties, "maxNoGridNodes", maxNoGridNodes);
  io::insert_property<SInd>(properties, "maxNoContainers", maxNoContainers);
  io::insert_property<std::function<void(Test3DGrid*)>>(properties, "meshGeneration", mesh_gen);
  io::insert_property<std::vector<Test3DGrid::Boundary>>(properties, "boundaries");

  return properties;
}

/// \test tests Grid constructor
TEST(node_container_test, constructor) {
  Test3DGrid grid(small_3D_grid());
  grid.read_mesh_generator();
  grid.generate_mesh();

  EXPECT_EQ(grid.nodes().size(),Ind(73));
}

Test3DGrid small3DGrid(small_3D_grid(),grid::initialize);



/// \test tests childs ranges
TEST(node_container_test, test_childs_range) {

  // test childs of root node
  auto expected_childs_0 = [](const SInd pos) {
    std::array<SInd,8> childs_0 = {{1,2,3,4,5,6,7,8}};
    return childs_0[pos];
  };
  auto childs = small3DGrid.nodes().childs(0);
  int posCounter = 0;
  for(auto child : childs) {
    EXPECT_EQ(child,expected_childs_0(posCounter));
    posCounter++;
  }

  // test childs of child 1
  auto expected_childs_1 = [](const SInd pos) {
    std::array<SInd,8> childs_1 = {{9,10,11,12,13,14,15,16}};
    return childs_1[pos];
  };
  childs = small3DGrid.nodes().childs(1);
  posCounter = 0;
  for(auto child : childs) {
    EXPECT_EQ(child,expected_childs_1(posCounter));
    posCounter++;
  }
}

/// \test tests child position in parent
TEST(node_container_test, test_child_position_in_parent) {

  // test for the root node
  auto expected_positions_0 = [](const SInt cId) {
    std::array<SInd,8> childs_1 = {{0,1,2,3,4,5,6,7}};
    return childs_1[cId - 1];
  };

  for(auto child : small3DGrid.nodes().childs(0)) {
    EXPECT_EQ(small3DGrid.nodes().position_in_parent(child),
              expected_positions_0(child));
  }

  // if the child has no parent you should expect a death in debug mode
  #ifndef NDEBUG
  EXPECT_DEATH(small3DGrid.nodes().position_in_parent(0),"[\\S\\s]+");
  #endif
}

/// \test leaf cell check
TEST(node_container_test, test_is_leaf_cell) {
  auto expected_leaf = [](const SInt nId) {
    return nId < 9 ? false : true;
  };
  for(auto nId : small3DGrid.nodes().nodes()) {
    EXPECT_EQ(small3DGrid.nodes().is_leaf(nId),expected_leaf(nId));
  }
}

/// \test test computation of samelvl nghbrs
TEST(node_container_test, test_samelvl_nghbrs) {

  /// check # of nghbrs returned is correct
  auto rootNghbrs = small3DGrid.nodes().all_samelvl_nghbrs(0);
  EXPECT_EQ(boost::distance(rootNghbrs), decltype(boost::distance(rootNghbrs))(2 * 3 /* 2 * nd*/));

  auto eId = invalid_value<Ind>();


  // checks nghbrs
  auto check_nghbrs = [](const Ind nId, const std::array<Ind,6>& nghbrsNid) {
    auto nghbrPos = small3DGrid.nodes().nghbr_pos();
    auto nghbrs = small3DGrid.nodes().all_samelvl_nghbrs(nId);
    for(auto pos : nghbrPos) {
      EXPECT_EQ(nghbrsNid[pos],nghbrs[pos]);
    }
  };

  auto nodeNghbrs = std::vector<std::array<Ind,6>>({
      {{ eId, eId, eId, eId, eId, eId }}, // nghbrs0
      {{ eId, 2  , eId, 3  , eId, 5   }}, // nghbrs1
      {{ 1  , eId, eId, 4  , eId, 6   }}, // nghbrs2
      {{ eId, 4  , 1  , eId, eId, 7   }}, // nghbrs3
      {{ 3  , eId, 2  , eId, eId, 8   }}, // nghbrs4
      {{ eId, 6  , eId, 7  , 1  , eId }}, // nghbrs5
      {{ 5  , eId, eId, 8  , 2  , eId }}, // nghbrs6
      {{ eId, 8  , 5  , eId, 3  , eId }}, // nghbrs7
      {{ 7  , eId, 6  , eId, 4  , eId }}, // nghbrs8
      {{ eId, 10 , eId, 11 , eId, 13  }}, // nghbrs9
      {{ 9  , 17 , eId, 12 , eId, 14  }}, // nghbrs10
      {{ eId,  12,   9, 25 , eId, 15  }} // nghbrs11
      // {{ 9  , eId, eId, 12 , eId, 14  }}, // nghbrs10
      // {{ 9  , eId, eId, 12 , eId, 14  }}, // nghbrs10
      // {{ 9  , eId, eId, 12 , eId, 14  }}, // nghbrs10
  });

 Ind nId = 0;
 for(auto nghbrs : nodeNghbrs) {
   check_nghbrs(nId,nghbrs);
   nId++;
 }

 // some nghbrs of some more complicated cells at the next level
 // auto nghbrs9  = std::array<Ind,6>{{ eId, 10 , eId, 11 , eId, 13  }};
 // auto nghbrs10 = std::array<Ind,6>{{ 9  , eId, eId, 4  , eId, 6   }};
 // auto nghbrs11 = std::array<Ind,6>{{ eId, 4  , 1  , eId, eId, 7   }};
 // auto nghbrs12 = std::array<Ind,6>{{ 3  , eId, 2  , eId, eId, 8   }};
 // node X is a corner of the domain
 // node Y cell whose parent and those of their neighbors are different
 // node Z cell whose parent and those of their neighbors are same

 // lets create a new grid and refine isotropically for one level
 // refine a single cell (creates a refinement jump)
 // test all neighbors

 // lets refine the other cells and refine again:
 // - a corner cell
 // - a cell in the middle (whose neighbors have different parents)
 // - test almost all neighbors
}

/// \test tests computation of cell neighbors including 2:1 jumps
TEST(node_container_test, test_nghbrs) {

}

/// \test tests nghbr ranges
TEST(node_container_test, test_nghbrs_ranges) {

}


////////////////////////////////////////////////////////////////////////////////
#undef ENABLE_DBG_
////////////////////////////////////////////////////////////////////////////////
