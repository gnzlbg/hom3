#ifndef HOM3_QUADRATURE_QUADRATURE_HPP_
#define HOM3_QUADRATURE_QUADRATURE_HPP_
////////////////////////////////////////////////////////////////////////////////
// Include:
#include "../globals.hpp"
#include "../grid/grid.hpp"
#include "../geometry/cell_cartesian.hpp"
////////////////////////////////////////////////////////////////////////////////

namespace quadrature {

template<SInd nd> struct GaussPoints {
  static const SInd size = Grid<nd>::no_edge_vertices();
  using type = std::array<NumA<nd>,size>;
  using dim = typename Grid<nd>::dim;

  type operator()() {
    return gauss_points_(dim());
  }

  type gauss_points_(grid::two_dim) {
    const NumA<nd> gp1 = { -1/std::sqrt(3), -1/std::sqrt(3) };
    const NumA<nd> gp2 = {  1/std::sqrt(3), -1/std::sqrt(3) };
    const NumA<nd> gp3 = { -1/std::sqrt(3),  1/std::sqrt(3) };
    const NumA<nd> gp4 = {  1/std::sqrt(3),  1/std::sqrt(3) };
    return {{gp1,gp2,gp3,gp4}};
  }

  type gauss_points_(grid::three_dim) {
    const NumA<nd> gp1 = { -1/std::sqrt(3), -1/std::sqrt(3), -1/std::sqrt(3) };
    const NumA<nd> gp2 = {  1/std::sqrt(3), -1/std::sqrt(3), -1/std::sqrt(3) };
    const NumA<nd> gp3 = { -1/std::sqrt(3),  1/std::sqrt(3), -1/std::sqrt(3) };
    const NumA<nd> gp4 = {  1/std::sqrt(3),  1/std::sqrt(3), -1/std::sqrt(3) };
    const NumA<nd> gp5 = { -1/std::sqrt(3), -1/std::sqrt(3),  1/std::sqrt(3) };
    const NumA<nd> gp6 = {  1/std::sqrt(3), -1/std::sqrt(3),  1/std::sqrt(3) };
    const NumA<nd> gp7 = { -1/std::sqrt(3),  1/std::sqrt(3),  1/std::sqrt(3) };
    const NumA<nd> gp8 = {  1/std::sqrt(3),  1/std::sqrt(3),  1/std::sqrt(3) };
    return {{gp1,gp2,gp3,gp4,gp5,gp6,gp7,gp8}};
  }
};

template<SInd nd, class Vector>
auto gauss_points_x(Vector&& x, const Num dx)
    -> decltype(GaussPoints<nd>()()){
  const Num dx2 = 0.5 * dx;
  auto gps_ref = GaussPoints<nd>()();
  for(auto& gp_ref : gps_ref) {
    gp_ref = x + dx2 * gp_ref;
  }
  return gps_ref;
};


template<SInd nd, class Vector, class Functor> auto integrate(Functor&& f, Vector&& xc, const Num dx)
    -> decltype(f(xc)) {
  auto x_gps = gauss_points_x<nd>(std::forward<Vector>(xc),dx);
  /// \todo stop hating c++ and use something better
  typename std::remove_reference<decltype(f(xc))>::type v_gp;
  std::array<decltype(v_gp),Grid<nd>::no_edge_vertices()> v_gps;
  for(SInd gp_i = 0; gp_i < Grid<nd>::no_edge_vertices(); ++gp_i) {
    v_gps[gp_i] = f(x_gps[gp_i]);
  }
  typename std::remove_reference<decltype(f(xc))>::type result = math::zero(decltype(f(xc))());
  for(SInd gp_i = 0; gp_i < Grid<nd>::no_edge_vertices(); ++gp_i) {
    result = result + v_gps[gp_i];
  }
  result = result * geometry::cell::cartesian::volume<nd>(dx) / math::ct::ipow(2,nd);
  return result;
}

} // quadrature namespace

////////////////////////////////////////////////////////////////////////////////
#endif
