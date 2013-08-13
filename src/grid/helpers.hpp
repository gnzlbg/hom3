#ifndef GRID_HELPERS_HPP_
#define GRID_HELPERS_HPP_
////////////////////////////////////////////////////////////////////////////////
/// Includes:
#include "../globals.hpp"
#include "grid.hpp"
#include "../geometry/geometry.hpp"
////////////////////////////////////////////////////////////////////////////////
namespace hom3 { namespace grid {
////////////////////////////////////////////////////////////////////////////////

namespace helpers {

/// \brief Easily builds a nd-hypercube grid
namespace cube {

/// \brief #of leaf nodes at \p level
template<SInd nd> static constexpr Ind no_leaf_nodes(const SInd level) {
  return math::ct::ipow(math::ct::ipow(2,nd),level);
}

/// \brief total #of nodes from \p level [0..level]
template<SInd nd> static constexpr Ind no_nodes(const SInd level) {
  return level == 0 ? 1 : no_leaf_nodes<nd>(level) + no_nodes<nd>(level-1);
}

/// \brief #of leaf nodes at each nd-hypercube face
template<SInd nd> static constexpr Ind no_nodes_per_face(const SInd level) {
  return math::ct::ipow(2, level*(nd - 1));
}

/// \brief #of nd-hypercube faces
template<SInd nd> static constexpr Ind no_faces() { return 2 * nd; }

/// \brief #of cells in solver container (including ghost cells)
/// N = #of_leafs + #of_faces * #of_cells_per_face
template<SInd nd> static constexpr Ind no_solver_cells_with_gc(const SInd level) {
  return no_leaf_nodes<nd>(level) // leaf cells
      + no_faces<nd>() * no_nodes_per_face<nd>(level); // ghost cells
}


namespace edge {

struct neg {}; struct pos {}; // negative/positive boundary
template<SInd nd> struct x {}; // perpendicular to the nd-axis
template<SInd nd> struct xneg : x<nd>, neg {}; // tag refinements
template<SInd nd> struct xpos : x<nd>, pos {};

/// Signs of the boundary surfaces normal vectors
constexpr Num n(neg) { return +1; }
constexpr Num n(pos) { return -1; }

// build vector with diagonal d along the axis x<i> and offset value off note
///
// clang bug: "= traits::dummy" is a hack, replacing it with pack expansion
// "..."  should work (in gcc it works, and standarese says it should work, and
// there is a bug filled in clang:
// http://llvm.org/bugs/show_bug.cgi?id=11723). The problem is that if you tell
// clang to expand an unused parameter pack it will just ignore it breaking
// SFNIAE which results in an ambiguity and the bug was filled in 2012 and noone
// is assigned yet and this month clang announced full c++11 support and that
// makes me angry and I just hate them.
template<SInd nd, traits::EnableIf<traits::equal<SInd,nd,2>> = traits::dummy >
NumA<nd> vector(x<0>,Num d,Num off){ return NumA<nd>{d,off}; }
template<SInd nd, traits::EnableIf<traits::equal<SInd,nd,3>> = traits::dummy >
NumA<nd> vector(x<0>,Num d,Num off){ return NumA<nd>{d,off,off}; }
template<SInd nd, traits::EnableIf<traits::equal<SInd,nd,2>> = traits::dummy >
NumA<nd> vector(x<1>,Num d,Num off){ return NumA<nd>{off,d}; }
template<SInd nd, traits::EnableIf<traits::equal<SInd,nd,3>> = traits::dummy >
NumA<nd> vector(x<1>,Num d,Num off){ return NumA<nd>{off,d,off}; }
template<SInd nd, traits::EnableIf<traits::equal<SInd,nd,3>> = traits::dummy >
NumA<nd> vector(x<2>,Num d,Num off) {return {off,off,d};}

// build surface center point, and surface normal
template<SInd nd, class edg> NumA<nd> center(edg, Num d, Num off) {
 return vector<nd>(edg(),d,off);
}
template<SInd nd, class edg> NumA<nd> normal(edg) {
  return vector<nd>(edg(),n(edg()),0.0);
}

template<SInd nd, class edg> auto make_edge(edg, Num d, Num off)
->decltype(geometry::make_geometry<
           geometry::implicit::Edge<nd>
           >(center<nd>(edg(),d,off), normal<nd>(edg()))) {
  return geometry::make_geometry<
    geometry::implicit::Edge<nd>
    >(center<nd>(edg(),d,off), normal<nd>(edg()));
}

} // edge namespace

using edge::make_edge;

template<SInd nd,class position>
boundary::Interface<nd> make_boundary
(position, Num d, Num off, std::string name, SolverIdx solverIdx,
 std::vector<boundary::ConditionHandler> bcs) {
  return {name,make_edge<nd>(position(),d,off),solverIdx,std::move(bcs)};
}

template<SInd nd> struct make_conditions {
  template<class T, traits::EnableIf<traits::equal<SInd,nd,2,T>> = traits::dummy>
  std::tuple<T,T,T,T> operator()(T t) {
    return std::make_tuple(t,t,t,t);
  }
  template<class T, traits::EnableIf<traits::equal<SInd,nd,3,T>> = traits::dummy>
  std::tuple<T,T,T,T,T,T> operator()(T t) {
    return std::make_tuple(t,t,t,t,t,t);
  }
};

template<SInd nd> struct make_boundaries {
  using Boundary = typename Grid<nd>::Boundary;
  using Boundaries = typename Grid<nd>::Boundaries;
  using Condition = boundary::ConditionHandler;

  struct BoundaryPosition { Num x_max, x_min, off; };
  auto boundary_position(RootCell<nd> rootCell) -> BoundaryPosition {
    auto length = rootCell.length;
    BoundaryPosition pos;
    auto eps = length * 1e-15; // std::numeric_limits<Num>::min();
    pos.x_max = rootCell.x_max()(0) - eps;
    pos.x_min = rootCell.x_min()(0) + eps;
    pos.off = pos.x_min + 0.5*length;
    return pos;
  }

  template<class T>
  Boundaries two_dim(SolverIdx solverId, RootCell<nd> rootCell, T&& bcs) {
    auto p = boundary_position(rootCell);

    Boundaries boundaries;
    boundaries.push_back(make_boundary<nd>(edge::xneg<0>(),p.x_min,p.off,"left"  ,{solverId},{std::get<0>(bcs)}));
    boundaries.push_back(make_boundary<nd>(edge::xpos<0>(),p.x_max,p.off,"right" ,{solverId},{std::get<1>(bcs)}));
    boundaries.push_back(make_boundary<nd>(edge::xneg<1>(),p.x_min,p.off,"bottom",{solverId},{std::get<2>(bcs)}));
    boundaries.push_back(make_boundary<nd>(edge::xpos<1>(),p.x_max,p.off,"top"   ,{solverId},{std::get<3>(bcs)}));
    return boundaries;
  }

  template<class T, traits::EnableIf<traits::equal<SInd,nd,2,T>> = traits::dummy>
  Boundaries operator()(SolverIdx solverId, RootCell<nd> rootCell, T&& bcs) {
    return two_dim(solverId,rootCell,std::forward<T>(bcs));
  }

  template<class T, traits::EnableIf<traits::equal<SInd,nd,3,T>> = traits::dummy>
  Boundaries operator()(SolverIdx solverId, RootCell<nd> rootCell, T&& bcs) {
    auto p = boundary_position(rootCell);
    auto boundaries = two_dim(solverId, rootCell, std::forward<T>(bcs));
    boundaries.push_back(make_boundary<nd>(edge::xneg<2>(),p.x_min,p.off,"front" ,{solverId},{std::get<4>(bcs)}));
    boundaries.push_back(make_boundary<nd>(edge::xpos<2>(),p.x_max,p.off,"back " ,{solverId},{std::get<5>(bcs)}));
    return boundaries;
  }
};


template<SInd nd, class MeshGen = grid::generation::MinLevel>
io::Properties properties
(const RootCell<nd> rootCell, const SInd minRefLvl, const SInd maxNoContainers = 1) {

  const Ind maxNoGridNodes = no_nodes<nd>(minRefLvl);

  io::Properties meshGeneration;
  io::insert_property<Ind>(meshGeneration,"level",minRefLvl);
  MeshGen mesh_gen(meshGeneration);

  io::Properties properties;
  io::insert_property<grid::RootCell<nd>>(properties,"rootCell",rootCell);
  io::insert_property<Ind>(properties,"maxNoGridNodes",maxNoGridNodes);
  io::insert_property<std::function<void(Grid<nd>&)>>(properties,"meshGeneration", mesh_gen);
  io::insert_property<SInd>(properties,"maxNoContainers",maxNoContainers);
  return properties;
}

} // cube namespace

} // helpers namespace

////////////////////////////////////////////////////////////////////////////////
}} // hom3::grid namespace
////////////////////////////////////////////////////////////////////////////////
#endif
