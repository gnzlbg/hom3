/// Includes:
#include "../../../grid/grid.hpp"
/// External Includes:
#include "gtest/gtest.h"
/// Options:
#define ENABLE_DBG_ 0
#include "../../../misc/dbg.hpp"
////////////////////////////////////////////////////////////////////////////////

template<SInd nd> io::Properties small_grid(const SInd minRefLevel) {
  using Boundaries = typename Grid<nd>::Boundaries;
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
void consistency_nghbr_check(const Grid<nd>& grid) {
  for(const auto nId : grid.nodes().nodes()) {
    for(const auto nghbrPos : grid.nodes().nghbr_pos()) {
      const auto nghbrId = grid.nodes().find_samelvl_nghbr(nId,nghbrPos);
      if(is_valid(nghbrId)) {
        const auto oppositeNghbrPos = grid.nodes().opposite_nghbr_position(nghbrPos);
        EXPECT_EQ(nId,grid.nodes().find_samelvl_nghbr(nghbrId,oppositeNghbrPos));
      }
    }
  }
}

/// \test tests Grid constructor
TEST(hierarchical_container_test, constructor) {
  Grid<3> grid(small_grid<3>(2));
  grid.read_mesh_generator();
  grid.generate_mesh();

  EXPECT_EQ(grid.nodes().size(),Ind(73));
}

Grid<3> small3DGrid(small_grid<3>(2),grid::initialize);

TEST(hierarchical_container_test, test_write_grid_domain_3d) {
  write_domain("grid_3D",&small3DGrid);
}

/// \test tests childs ranges
TEST(hierarchical_container_test, test_childs_range) {

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
TEST(hierarchical_container_test, test_child_position_in_parent) {

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
TEST(hierarchical_container_test, test_is_leaf_cell) {
  auto expected_leaf = [](const SInt nId) {
    return nId < 9 ? false : true;
  };
  for(auto nId : small3DGrid.nodes().nodes()) {
    EXPECT_EQ(small3DGrid.nodes().is_leaf(nId),expected_leaf(nId));
  }
}

/// \test test computation of samelvl nghbrs
TEST(hierarchical_container_test, test_3D_strict_samelvl_nghbrs) {

  /// check # of nghbrs returned is correct
  auto rootNghbrs = small3DGrid.nodes().all_samelvl_nghbrs<strict>(0);
  EXPECT_EQ(boost::distance(rootNghbrs), decltype(boost::distance(rootNghbrs))(2 * 3 /* 2 * nd*/));

  consistency_nghbr_check(small3DGrid);

  auto eId = invalid_value<Ind>();


  // checks nghbrs
  auto check_nghbrs = [](const Ind nId, const std::array<Ind,6>& nghbrsNid) {
    auto nghbrPos = small3DGrid.nodes().nghbr_pos();
    auto nghbrs = small3DGrid.nodes().all_samelvl_nghbrs<strict>(nId);
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

/// \test test computation of samelvl nghbrs
TEST(hierarchical_container_test, test_lazy_samelvl_nghbrs) {
}

/// \test tests computation of cell neighbors including 2:1 jumps
TEST(hierarchical_container_test, test_nghbrs) {

}

/// \test tests nghbr ranges
TEST(hierarchical_container_test, test_nghbrs_ranges) {

}


Grid<2> small2DGrid(small_grid<2>(3),grid::initialize);

TEST(hierarchical_container_test, test_write_grid_domain_2d) {
  write_domain("grid_2D",&small2DGrid);
}

/// \test test computation of samelvl nghbrs
TEST(hierarchical_container_test, test_2D_strict_samelvl_nghbrs) {

  /// check # of nghbrs returned is correct
  auto rootNghbrs = small2DGrid.nodes().all_samelvl_nghbrs<strict>(0);
  EXPECT_EQ(boost::distance(rootNghbrs), decltype(boost::distance(rootNghbrs))(2 * 2 /* 2 * nd*/));

  auto eId = invalid_value<Ind>();

  consistency_nghbr_check(small2DGrid);

  // checks nghbrs
  auto check_nghbrs = [](const Ind nId, const std::array<Ind,4>& shouldNghbrId) {
    auto nghbrPos = small2DGrid.nodes().nghbr_pos();
    auto isNghbrId = small2DGrid.nodes().all_samelvl_nghbrs<strict>(nId);
    //DBGV_ON((nId)(isNghbrId[0])(shouldNghbrId[0])(isNghbrId[0])(shouldNghbrId[1])(isNghbrId[1])(shouldNghbrId[2])(isNghbrId[2])(shouldNghbrId[3])(isNghbrId[3]));
    for(auto pos : nghbrPos) {
      EXPECT_EQ(shouldNghbrId[pos],isNghbrId[pos]);
    }
  };


  // expected output
  auto nodeNghbrs = std::vector<nghbrIds<2>>({
      // nId, [nghbr0...nghbr3]
      {  0, {{ eId, eId, eId, eId }}},
      {  1, {{ eId, 2  , eId, 3   }}},
      {  2, {{ 1  , eId, eId, 4   }}},
      {  3, {{ eId, 4  , 1  , eId }}},
      {  4, {{ 3  , eId, 2  , eId }}},
      {  5, {{ eId, 6  , eId, 7   }}},
      {  6, {{ 5  , 9  , eId, 8   }}},
      {  7, {{ eId, 8  , 5  , 13  }}},
      {  8, {{ 7  , 11 , 6  , 14  }}},
      {  9, {{ 6  , 10 , eId, 11  }}},
      { 10, {{ 9  , eId, eId, 12  }}},
      { 11, {{ 8  ,  12,   9, 17  }}},
      { 12, {{ 11 , eId, 10 , 18  }}},
      { 13, {{ eId, 14 , 7  , 15  }}},
      { 14, {{ 13 , 17 , 8  , 16  }}},
      { 15, {{ eId, 16 , 13 , eId }}},
      { 16, {{ 15 , 19 , 14 , eId }}},
      { 17, {{ 14 , 18 , 11 , 19  }}},
      { 18, {{ 17 , eId, 12 , 20  }}},
      { 19, {{ 16 , 20 , 17 , eId }}},
      { 20, {{ 19 , eId, 18 , eId }}},
      { 21, {{ eId, 22 , eId, 23  }}},
      { 24, {{ 23 , 27 , 22 , 30  }}},
      { 27, {{ 24 , 28 , 25 , 33  }}},
      { 28, {{ 27 , 39 , 26 , 34  }}},
      { 30, {{ 29 , 33 , 24 , 32  }}},
      { 32, {{ 31 , 35 , 30 , 54  }}},
      { 33, {{ 30 , 34 , 27 , 35  }}},
      { 34, {{ 33 , 45 , 28 , 36  }}},
      { 35, {{ 32 , 36 , 33 , 57  }}},
      { 36, {{ 35 , 47 , 34 , 58  }}},
      { 39, {{ 28 , 40 , 37 , 45  }}},
      { 40, {{ 39 , 43 , 38 , 46  }}},
      { 42, {{ 41 , eId, eId, 44  }}},
      { 43, {{ 40 , 44 , 41 , 49  }}},
      { 45, {{ 34 , 46 , 39 , 47  }}},
      { 46, {{ 45 , 49 , 40 , 48  }}},
      { 47, {{ 36 , 48 , 45 , 69  }}},
      { 48, {{ 47 , 51 , 46 , 70  }}},
      { 49, {{ 46 , 50 , 43 , 51  }}},
      { 54, {{ 53 , 57 , 32 , 56  }}},
      { 57, {{ 54 , 58 , 35 , 59  }}},
      { 58, {{ 57 , 69 , 36 , 60  }}},
      { 69, {{ 58 , 70 , 47 , 71  }}},
      { 70, {{ 69 , 73 , 48 , 72  }}},
      { 73, {{ 70 , 74 , 51 , 75  }}},
      { 56, {{ 55 , 59 , 54 , 62  }}},
      { 59, {{ 56 , 60 , 57 , 65  }}},
      { 60, {{ 59 , 71 , 58 , 66  }}},
      { 71, {{ 60 , 72 , 69 , 77  }}},
      { 72, {{ 71 , 75 , 70 , 78  }}},
      { 73, {{ 70 , 74 , 51 , 75  }}},
      { 62, {{ 61 , 65 , 56 , 64  }}},
      { 65, {{ 62 , 66 , 59 , 67  }}},
      { 66, {{ 65 , 77 , 60 , 68  }}},
      { 77, {{ 66 , 78 , 71 , 79  }}},
      { 78, {{ 77 , 81 , 72 , 80  }}},
      { 81, {{ 78 , 82 , 75 , 83  }}},
      { 63, {{ eId, 64 , 61 , eId }}},
      { 84, {{ 83 , eId, 82 , eId }}}
  });

  for(auto nghbrs : nodeNghbrs) {
   check_nghbrs(nghbrs.nId,nghbrs.nghbrs);
 }
}

////////////////////////////////////////////////////////////////////////////////
#undef ENABLE_DBG_
////////////////////////////////////////////////////////////////////////////////
