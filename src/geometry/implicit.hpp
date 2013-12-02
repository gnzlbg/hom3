#ifndef HOM3_GEOMETRY_IMPLICIT_HPP_
#define HOM3_GEOMETRY_IMPLICIT_HPP_
////////////////////////////////////////////////////////////////////////////////
#include <type_traits>
#include <algorithm>
/// \file \brief Contains functions for describing implicit geometries
///
/// \todo rename this file to implicit.hpp
/// \todo: does using homogenized coordinates makes sense? (don't think so)
/// \todo: refactor rotations as an adaptor
/// \todo: does adaptor::scale makes sense? (don't think so)
/// \todo: think about performance cost of adaptors (kind of done, combine takes
/// a functor now). Still, accessing geometries by pointer might not be the
/// best.
////////////////////////////////////////////////////////////////////////////////
namespace hom3 {
////////////////////////////////////////////////////////////////////////////////

/// \brief Geometry-modelling utilities
namespace geometry {

/// \brief Implements the geometry concept
template<SInd nd> struct Interface {
  template<class G> explicit Interface(G&& g) noexcept
    : signed_distance([&](const NumA<nd>& x) { return g(x); }) {}
  const std::function<Num(const NumA<nd>&)> signed_distance;
};

template<class Geometry, class... Args>
std::shared_ptr<Geometry> make_geometry(Args&&... args)
{ return std::make_shared<Geometry>(std::forward<Args>(args)...); }

/// \brief Implicit signed-distance functions
///
/// Convention: g(x) > 0 outside, g(x) = 0 at surface, g(x) < 0 inside
namespace implicit {

/// \brief Signed-distance to sphere of radius r centered at xc = (x,y,z)^T.
template<SInd nd> struct Sphere {
  static const SInd no_dims = nd;
  const NumA<nd> xc;
  const Num r;
  Sphere(const NumA<nd>& center, Num radius) noexcept : xc(center), r(radius) {}
  Num operator()(const NumA<nd>& x) const noexcept
  { return (x - xc).norm() - r; }
};

namespace moving {

template<SInd nd> struct Sphere {
  static const SInd no_dims = nd;
  std::function<NumA<nd>()> xc;
  const Num r;
  Sphere(std::function<NumA<nd>()> center, Num radius) noexcept : xc(center), r(radius) {}
  Num operator()(const NumA<nd>& x) const noexcept
  { return (x - xc()).norm() - r; }
};

}  // moving

/// \brief Signed-distance to hexaedron of lengths L=(lx,ly,lz)^T,
/// centered at xc = (x,y,z)^T, and rotated by phi = (phi_x,phi_y,phi_z)^T.
///
/// \todo Implement rotations
template<SInd nd> struct Square {
  static const SInd no_dims = nd;
  const NumA<nd> xc;
  const NumA<nd> l_2;
  Square(const NumA<nd>& xc_, NumA<nd> l_) noexcept : xc(xc_), l_2(0.5*l_) {}
  Num operator()(const NumA<nd>& x) const noexcept {
    // compute signed distances per dir
    NumA<nd> distance = xc - l_2 - x;  // x(d) < xc(d);
    for (SInd d = 0; d < nd; ++d) {
      if (x(d) > xc(d)) {
        distance(d) = x(d) - (xc(d) + l_2(d));  // pos if outside, neg if inside
      }
    }

    // #of positive values determines the quadrant type
    const SInd noPosVals = (distance.array() > 0).count();

    if (noPosVals == nd) {  // corner quadrants
      return distance.norm();
    } else if (noPosVals < 2) {  // direct neighbors and inside
      return distance.maxCoeff();
    } else if (nd == 3) {  // shorter distanceance is to edges
      if (distance(0) < 0) {  // parallel x
        return std::sqrt(std::pow(distance(1), 2) + std::pow(distance(2), 2));
      } else if (distance(1) < 0) {  // parallel y
        return std::sqrt(std::pow(distance(0), 2) + std::pow(distance(2), 2));
      } else if (distance(2) < 0) {  // parallel z
        return std::sqrt(std::pow(distance(0), 2) + std::pow(distance(1), 2));
      }
    }

    TERMINATE("This function should exit before");
  }
};

/// \brief Signed-disntance to an edge. Edges are oriented:
/// - points in 1D,
/// - lines in 2D, and
/// - planes in 3D.
template<SInd nd> struct Edge {
  static const SInd no_dims = nd;
  const NumA<nd> point;
  /// Normalized normal vector
  const NumA<nd> normal;
  /// Direction (only used for 2D)
  const NumA<nd> dir;

  Edge(const NumA<nd> p, const NumA<nd> n) noexcept : point(p), normal(n) {}

  inline Num operator()(const NumA<nd>& x) const noexcept
  { return normal.dot(x - point); }
};

////////////////////////////////////////////////////////////////////////////////
/// \brief Implicit function adaptors
namespace adaptors {
////////////////////////////////////////////////////////////////////////////////

/// \brief Inverts the sign of an implicit function
template<class T> struct Invert : T {
  template<class...Ts> Invert(Ts&&... ts) noexcept
    : T(std::forward<Ts>(ts)...) {}

  Num operator()(const NumA<T::no_dims>& x) const noexcept
  { return -1. * T::operator()(x); }
};

template<class T> struct Inverted {
  explicit Inverted(std::shared_ptr<T> t) noexcept : t_(t) {}

  Num operator()(const NumA<T::no_dims>& x) const noexcept
  { return -1 * t_->operator()(x); }

 private:
  std::shared_ptr<T> t_;
};
template<class T> Inverted<T> invert2(std::shared_ptr<T> t) noexcept
{ return Inverted<T>(t); }
template<class T>
std::shared_ptr<Inverted<T>> invert(std::shared_ptr<T> t) noexcept
{ return make_geometry<Inverted<T>>(t); }

template<class T, class U, class C> struct Combine {
  Combine(std::shared_ptr<T> t, std::shared_ptr<U> u) noexcept : t_(t), u_(u) {}
  Num operator()(const NumA<T::no_dims>& x) const noexcept {
    static_assert(T::no_dims == U::no_dims,
                  "Geometry binary union: dimension mismatch!");
    const Num d1 = t_->operator()(x);
    const Num d2 = u_->operator()(x);
    return C()(d1, d2);  // maybe its better to create a data member C c_?
  }
 private:
  std::shared_ptr<T> t_;
  std::shared_ptr<U> u_;
};

struct UnionF {
  Num operator()(const Num d1, const Num d2) const noexcept
  { return std::min(d1, d2); }
};
struct IntersectionF {
  Num operator()(const Num d1, const Num d2) const noexcept
  { return std::max(d1, d2); }
};
struct DifferenceF {
  Num operator()(const Num d1, const Num d2) const noexcept
  { return std::max(d1, -d2); }
};

template<class T, class U> using Union        = Combine<T, U, UnionF>;
template<class T, class U> using Intersection = Combine<T, U, IntersectionF>;
template<class T, class U> using Difference   = Combine<T, U, DifferenceF>;

template<class T, class U>
Union<T, U> make_union(std::shared_ptr<T> t, std::shared_ptr<U> u) noexcept
{ return Union<T, U>(t, u); }
template<class T, class U>
Intersection<T, U> make_intersection(std::shared_ptr<T> t, std::shared_ptr<U> u)
{ return Intersection<T, U>(t, u); }
template<class T, class U>
Difference<T, U> make_difference(std::shared_ptr<T> t, std::shared_ptr<U> u)
{ return Difference<T, U>(t, u); }

/// Translate: should enrich interface
/// Rotate : should enrich interface

}  // namespace adaptors
}  // namespace implicit

/// \brief Makes a cut-off cube
template<SInd nd> auto make_cube
(const NumA<nd> x_center, const NumA<nd> dimensions0, const Num cell_length,
const Num tol = 0.6) noexcept {
  const NumA<nd> dimensions1 = dimensions0.array() + (cell_length);
  const NumA<nd> dimensions2 = dimensions0.array() - tol * (cell_length);

  auto cube0
    = geometry::make_geometry<
        geometry::implicit::Square<nd>
      >(x_center, dimensions0);

  auto cube1
    = geometry::make_geometry<
        geometry::implicit::Square<nd>
      >(x_center, dimensions1);

  auto cube2
    = geometry::implicit::adaptors::invert
      (geometry::make_geometry<
         geometry::implicit::Square<nd>
       >(x_center, dimensions2));

  return std::make_tuple(cube0, cube1, cube2);
}

////////////////////////////////////////////////////////////////////////////////
}  // namespace geometry
}  // namespace hom3
////////////////////////////////////////////////////////////////////////////////
#endif
