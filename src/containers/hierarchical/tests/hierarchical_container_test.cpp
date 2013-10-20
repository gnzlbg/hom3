/// Includes:
#include "grid/grid.hpp"
#include "grid/helpers.hpp"
/// External Includes:
#include "misc/test.hpp"
/// Options:
#define ENABLE_DBG_ 0
#include "misc/dbg.hpp"
////////////////////////////////////////////////////////////////////////////////

using namespace hom3;

template<SInd nd> io::Properties small_grid(const SInd minRefLevel) {
  using Boundaries = typename grid::Grid<nd>::Boundaries;
  Boundaries boundaries;
  auto properties
      = grid::helpers::cube::properties<nd>
      ({NumA<nd>::Constant(0),NumA<nd>::Constant(1)}, minRefLevel);
  io::insert_property<decltype(boundaries)>(properties,"boundaries",boundaries);
  return properties;
}

template<SInd nd> struct nghbrIds {
  Ind nId;
  std::array<Ind,2*nd> nghbrs;
};

/// verify that if a node nId has another node nghbrId as nghbr in pos nghbrPos,
/// nghbrId must have nId as neighbor also but in position oppositeNghbrPos
template<SInd nd>
void consistency_nghbr_check(const grid::Grid<nd>& grid) {
  for(const auto nId : grid.nodes()) {
    for(const auto nghbrPos : grid.neighbor_positions()) {
      const auto nghbrId = grid.find_samelvl_neighbor(nId,nghbrPos);
      if(is_valid(nghbrId)) {
        const auto oppositeNghbrPos = grid.opposite_neighbor_position(nghbrPos);
        EXPECT_EQ(nId,grid.find_samelvl_neighbor(nghbrId,oppositeNghbrPos));
      }
    }
  }
}

/// \test tests Grid constructor
TEST(hierarchical_container_test, constructor) {
  grid::Grid<3> grid(small_grid<3>(2));
  grid.read_mesh_generator();
  grid.generate_mesh();

  EXPECT_EQ(grid.size(),Ind(73));
}

grid::Grid<3> small3DGrid(small_grid<3>(2),grid::initialize);

TEST(hierarchical_container_test, test_write_grid_domain_3d) {
  write_domain("grid_3D", small3DGrid);
}

/// \test tests childs ranges
TEST(hierarchical_container_test, test_childs_range) {

  // test childs of root node
  auto expected_childs_0 = [](const SInd pos) {
    std::array<SInd,8> childs_0 = {{1,2,3,4,5,6,7,8}};
    return childs_0[pos];
  };
  auto childs = small3DGrid.childs(NodeIdx{0});
  int posCounter = 0;
  for(auto child : childs) {
    EXPECT_EQ(expected_childs_0(posCounter),child());
    posCounter++;
  }

  // test childs of child 1
  auto expected_childs_1 = [](const SInd pos) {
    std::array<SInd,8> childs_1 = {{9,10,11,12,13,14,15,16}};
    return childs_1[pos];
  };
  childs = small3DGrid.childs(NodeIdx{1});
  posCounter = 0;
  for(auto child : childs) {
    EXPECT_EQ(expected_childs_1(posCounter),child());
    posCounter++;
  }
}

/// \test tests child position in parent
TEST(hierarchical_container_test, test_child_position_in_parent) {

  // test for the root node
  auto expected_positions_0 = [](const NodeIdx nIdx) {
    std::array<SInd,8> childs_1 = {{0,1,2,3,4,5,6,7}};
    return childs_1[(nIdx - NodeIdx{1})()];
  };

  for(auto child : small3DGrid.childs(NodeIdx{0})) {
    EXPECT_EQ(expected_positions_0(child),small3DGrid.position_in_parent(child));
  }

  // if the child has no parent you should expect a death in debug mode
  #ifndef NDEBUG
  EXPECT_DEATH(small3DGrid.position_in_parent(NodeIdx{0}),"[\\S\\s]+");
  #endif
}

/// \test leaf cell check
TEST(hierarchical_container_test, test_is_leaf_cell) {
  auto expected_leaf = [](const NodeIdx nIdx) {
    return nIdx < NodeIdx{9} ? false : true;
  };
  for(auto nIdx : small3DGrid.nodes()) {
    EXPECT_EQ(expected_leaf(nIdx),small3DGrid.is_leaf(nIdx));
  }
}

/// \test test computation of samelvl nghbrs
TEST(hierarchical_container_test, test_3D_strict_samelvl_nghbrs) {

  /// check # of nghbrs returned is correct
  auto rootNghbrs = small3DGrid.all_samelvl_neighbors<strict>(NodeIdx{0});
  EXPECT_EQ(decltype(boost::distance(rootNghbrs))(2 * 3 /* 2 * nd*/),boost::distance(rootNghbrs));

  consistency_nghbr_check(small3DGrid);

  auto eIdx = invalid<Ind>();


  // checks nghbrs
  auto check_nghbrs = [](const NodeIdx nId, const std::array<Ind,6>& nghbrsNid) {
    auto nghbrPos = small3DGrid.neighbor_positions();
    auto nghbrs = small3DGrid.all_samelvl_neighbors<strict>(nId);
    for(auto pos : nghbrPos) {
      EXPECT_EQ(NodeIdx{nghbrsNid[pos]},nghbrs[pos]);
    }
  };

  auto nodeNghbrs = std::vector<std::array<Ind,6>>({
      {{ eIdx, eIdx, eIdx, eIdx, eIdx, eIdx }}, // nghbrs0
      {{ eIdx, 2  , eIdx, 3  , eIdx, 5   }}, // nghbrs1
      {{ 1  , eIdx, eIdx, 4  , eIdx, 6   }}, // nghbrs2
      {{ eIdx, 4  , 1  , eIdx, eIdx, 7   }}, // nghbrs3
      {{ 3  , eIdx, 2  , eIdx, eIdx, 8   }}, // nghbrs4
      {{ eIdx, 6  , eIdx, 7  , 1  , eIdx }}, // nghbrs5
      {{ 5  , eIdx, eIdx, 8  , 2  , eIdx }}, // nghbrs6
      {{ eIdx, 8  , 5  , eIdx, 3  , eIdx }}, // nghbrs7
      {{ 7  , eIdx, 6  , eIdx, 4  , eIdx }}, // nghbrs8
      {{ eIdx, 10 , eIdx, 11 , eIdx, 13  }}, // nghbrs9
      {{ 9  , 17 , eIdx, 12 , eIdx, 14  }}, // nghbrs10
      {{ eIdx,  12,   9, 25 , eIdx, 15  }} // nghbrs11
      // {{ 9  , eIdx, eIdx, 12 , eIdx, 14  }}, // nghbrs10
      // {{ 9  , eIdx, eIdx, 12 , eIdx, 14  }}, // nghbrs10
      // {{ 9  , eIdx, eIdx, 12 , eIdx, 14  }}, // nghbrs10
  });

 auto nId = NodeIdx{0};
 for(auto nghbrs : nodeNghbrs) {
   check_nghbrs(nId,nghbrs);
   ++nId;
 }

 // some nghbrs of some more complicated cells at the next level
 // auto nghbrs9  = std::array<Ind,6>{{ eIdx, 10 , eIdx, 11 , eIdx, 13  }};
 // auto nghbrs10 = std::array<Ind,6>{{ 9  , eIdx, eIdx, 4  , eIdx, 6   }};
 // auto nghbrs11 = std::array<Ind,6>{{ eIdx, 4  , 1  , eIdx, eIdx, 7   }};
 // auto nghbrs12 = std::array<Ind,6>{{ 3  , eIdx, 2  , eIdx, eIdx, 8   }};
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

/// \test test computation of samelvl nghbrs
TEST(hierarchical_container_test, test_lazy_samelvl_nghbrs) {
}

/// \test tests computation of cell neighbors including 2:1 jumps
TEST(hierarchical_container_test, test_nghbrs) {

}

/// \test tests nghbr ranges
TEST(hierarchical_container_test, test_nghbrs_ranges) {

}


grid::Grid<2> small2DGrid(small_grid<2>(3), grid::initialize);

TEST(hierarchical_container_test, test_write_grid_domain_2d) {
  write_domain("grid_2D", small2DGrid);
}

/// \test test computation of samelvl nghbrs
TEST(hierarchical_container_test, test_2D_strict_samelvl_nghbrs) {

  /// check # of nghbrs returned is correct
  auto rootNghbrs = small2DGrid.all_samelvl_neighbors<strict>(NodeIdx{0});
  EXPECT_EQ(boost::distance(rootNghbrs), decltype(boost::distance(rootNghbrs))(2 * 2 /* 2 * nd*/));

  auto eIdx = invalid<Ind>();

  consistency_nghbr_check(small2DGrid);

  // checks nghbrs
  auto check_nghbrs = [](const NodeIdx nIdx, const std::array<Ind,4>& shouldNghbrId) {
    auto nghbrPos = small2DGrid.neighbor_positions();
    auto isNghbrId = small2DGrid.all_samelvl_neighbors<strict>(nIdx);
    //DBGV_ON((nId)(isNghbrId[0])(shouldNghbrId[0])(isNghbrId[0])(shouldNghbrId[1])(isNghbrId[1])(shouldNghbrId[2])(isNghbrId[2])(shouldNghbrId[3])(isNghbrId[3]));
    for(auto pos : nghbrPos) {
      EXPECT_EQ(NodeIdx{shouldNghbrId[pos]},isNghbrId[pos]);
    }
  };


  // expected output
  auto nodeNghbrs = std::vector<nghbrIds<2>>({
      // nId, [nghbr0...nghbr3]
      {  0, {{ eIdx, eIdx, eIdx, eIdx }}},
      {  1, {{ eIdx, 2   , eIdx, 3    }}},
      {  2, {{ 1   , eIdx, eIdx, 4    }}},
      {  3, {{ eIdx, 4   , 1   , eIdx }}},
      {  4, {{ 3   , eIdx, 2   , eIdx }}},
      {  5, {{ eIdx, 6   , eIdx, 7    }}},
      {  6, {{ 5   , 9   , eIdx, 8    }}},
      {  7, {{ eIdx, 8   , 5   , 13   }}},
      {  8, {{ 7   , 11  , 6   , 14   }}},
      {  9, {{ 6   , 10  , eIdx, 11   }}},
      { 10, {{ 9   , eIdx, eIdx, 12   }}},
      { 11, {{ 8   , 12  , 9   , 17   }}},
      { 12, {{ 11  , eIdx, 10  , 18   }}},
      { 13, {{ eIdx, 14  , 7   , 15   }}},
      { 14, {{ 13  , 17  , 8   , 16   }}},
      { 15, {{ eIdx, 16  , 13  , eIdx }}},
      { 16, {{ 15  , 19  , 14  , eIdx }}},
      { 17, {{ 14  , 18  , 11  , 19   }}},
      { 18, {{ 17  , eIdx, 12  , 20   }}},
      { 19, {{ 16  , 20  , 17  , eIdx }}},
      { 20, {{ 19  , eIdx, 18  , eIdx }}},
      { 21, {{ eIdx, 22  , eIdx, 23   }}},
      { 24, {{ 23  , 27  , 22  , 30   }}},
      { 27, {{ 24  , 28  , 25  , 33   }}},
      { 28, {{ 27  , 39  , 26  , 34   }}},
      { 30, {{ 29  , 33  , 24  , 32   }}},
      { 32, {{ 31  , 35  , 30  , 54   }}},
      { 33, {{ 30  , 34  , 27  , 35   }}},
      { 34, {{ 33  , 45  , 28  , 36   }}},
      { 35, {{ 32  , 36  , 33  , 57   }}},
      { 36, {{ 35  , 47  , 34  , 58   }}},
      { 39, {{ 28  , 40  , 37  , 45   }}},
      { 40, {{ 39  , 43  , 38  , 46   }}},
      { 42, {{ 41  , eIdx, eIdx, 44   }}},
      { 43, {{ 40  , 44  , 41  , 49   }}},
      { 45, {{ 34  , 46  , 39  , 47   }}},
      { 46, {{ 45  , 49  , 40  , 48   }}},
      { 47, {{ 36  , 48  , 45  , 69   }}},
      { 48, {{ 47  , 51  , 46  , 70   }}},
      { 49, {{ 46  , 50  , 43  , 51   }}},
      { 54, {{ 53  , 57  , 32  , 56   }}},
      { 57, {{ 54  , 58  , 35  , 59   }}},
      { 58, {{ 57  , 69  , 36  , 60   }}},
      { 69, {{ 58  , 70  , 47  , 71   }}},
      { 70, {{ 69  , 73  , 48  , 72   }}},
      { 73, {{ 70  , 74  , 51  , 75   }}},
      { 56, {{ 55  , 59  , 54  , 62   }}},
      { 59, {{ 56  , 60  , 57  , 65   }}},
      { 60, {{ 59  , 71  , 58  , 66   }}},
      { 71, {{ 60  , 72  , 69  , 77   }}},
      { 72, {{ 71  , 75  , 70  , 78   }}},
      { 73, {{ 70  , 74  , 51  , 75   }}},
      { 62, {{ 61  , 65  , 56  , 64   }}},
      { 65, {{ 62  , 66  , 59  , 67   }}},
      { 66, {{ 65  , 77  , 60  , 68   }}},
      { 77, {{ 66  , 78  , 71  , 79   }}},
      { 78, {{ 77  , 81  , 72  , 80   }}},
      { 81, {{ 78  , 82  , 75  , 83   }}},
      { 63, {{ eIdx, 64  , 61  , eIdx }}},
      { 84, {{ 83  , eIdx, 82  , eIdx }}}
  });

  for(auto nghbrs : nodeNghbrs) {
    check_nghbrs(NodeIdx{nghbrs.nId},nghbrs.nghbrs);
 }
}

////////////////////////////////////////////////////////////////////////////////
#undef ENABLE_DBG_
////////////////////////////////////////////////////////////////////////////////
