#ifndef GRID_HELPERS_HPP_
#define GRID_HELPERS_HPP_
////////////////////////////////////////////////////////////////////////////////
/// Includes:
#include <string>
#include <vector>
#include "globals.hpp"
#include "grid.hpp"
#include "geometry/geometry.hpp"
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
template<SInd nd, EnableIf<traits::equal<SInd,nd,2>> = traits::dummy >
NumA<nd> vector(x<0>,Num d,Num off){ return NumA<nd>{d,off}; }
template<SInd nd, EnableIf<traits::equal<SInd,nd,3>> = traits::dummy >
NumA<nd> vector(x<0>,Num d,Num off){ return NumA<nd>{d,off,off}; }
template<SInd nd, EnableIf<traits::equal<SInd,nd,2>> = traits::dummy >
NumA<nd> vector(x<1>,Num d,Num off){ return NumA<nd>{off,d}; }
template<SInd nd, EnableIf<traits::equal<SInd,nd,3>> = traits::dummy >
NumA<nd> vector(x<1>,Num d,Num off){ return NumA<nd>{off,d,off}; }
template<SInd nd, EnableIf<traits::equal<SInd,nd,3>> = traits::dummy >
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

/// \brief Makes a single boundary condition
template<SInd nd, class Solver, class position, class BC>
typename Solver::Boundary make_boundary
(position, Num d, Num off, std::string name, Solver& solver, BC&& bc) {
  return {name, make_edge<nd>(position(), d, off), solver, std::forward<BC>(bc)};
}

template<SInd nd, class T, EnableIf<traits::equal<SInd,nd,2,T>> = traits::dummy>
std::tuple<T,T,T,T> make_conditions(T t) {
  return std::make_tuple(t,t,t,t);
}

template<SInd nd, class T, EnableIf<traits::equal<SInd,nd,3,T>> = traits::dummy>
std::tuple<T,T,T,T,T,T> make_conditions(T t) {
  return std::make_tuple(t,t,t,t,t,t);
}

/// \brief Build boundaries for an unit cube
template<SInd nd> struct make_boundaries {
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

  template<class Solver, class BCs>
  typename Solver::Boundaries
  two_dim(Solver& solver, RootCell<nd> rootCell, BCs&& bcs) {
    auto p = boundary_position(rootCell);

    typename Solver::Boundaries boundaries;
    boundaries.emplace_back
    (make_boundary<nd>(edge::xneg<0>(), p.x_min, p.off, "left"  , solver, std::get<0>(bcs)));
    boundaries.emplace_back
    (make_boundary<nd>(edge::xpos<0>(), p.x_max, p.off, "right" , solver, std::get<1>(bcs)));
    boundaries.emplace_back
    (make_boundary<nd>(edge::xneg<1>(), p.x_min, p.off, "bottom", solver, std::get<2>(bcs)));
    boundaries.emplace_back
    (make_boundary<nd>(edge::xpos<1>(), p.x_max, p.off, "top"   , solver, std::get<3>(bcs)));
    return boundaries;
  }

  template<class Solver, class BCs,
           EnableIf<traits::equal<SInd, nd, 2, BCs>> = traits::dummy>
  typename Solver::Boundaries
  operator()(Solver& solver, RootCell<nd> rootCell, BCs&& bcs) {
    return two_dim(solver, rootCell, std::forward<BCs>(bcs));
  }

  template<class Solver, class BCs,
           EnableIf<traits::equal<SInd, nd, 3, BCs>> = traits::dummy>
  typename Solver::Boundaries
  operator()(Solver& solver, RootCell<nd> rootCell, BCs&& bcs) {
    auto p = boundary_position(rootCell);
    auto boundaries = two_dim(solver, rootCell, std::forward<BCs>(bcs));
    boundaries.emplace_back
    (make_boundary<nd>(edge::xneg<2>(), p.x_min, p.off, "front", solver, std::get<4>(bcs)));
    boundaries.emplace_back
    (make_boundary<nd>(edge::xpos<2>(), p.x_max, p.off, "back ", solver, std::get<5>(bcs)));
    return boundaries;
  }
};


template<SInd nd, class MeshGen = grid::generation::MinLevel>
io::Properties properties
(const RootCell<nd> rootCell, const SInd minRefLvl, const SInd maxNoGridSolvers = 1) {

  const Ind maxNoGridNodes = no_nodes<nd>(minRefLvl);

  io::Properties meshGeneration;
  io::insert_property<Ind>(meshGeneration,"level",minRefLvl);
  MeshGen mesh_gen(meshGeneration);

  io::Properties properties;
  io::insert_property<grid::RootCell<nd>>(properties,"rootCell",rootCell);
  io::insert_property<Ind>(properties,"maxNoGridNodes",maxNoGridNodes);
  io::insert_property<std::function<void(Grid<nd>&)>>(properties,"meshGeneration", mesh_gen);
  io::insert_property<SInd>(properties,"maxNoGridSolvers",maxNoGridSolvers);
  return properties;
}

} // cube namespace

} // helpers namespace

////////////////////////////////////////////////////////////////////////////////
}} // hom3::grid namespace
////////////////////////////////////////////////////////////////////////////////
#endif
