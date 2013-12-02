#ifndef HOM3_GEOMETRY_ALGORITHMS_HPP_
#define HOM3_GEOMETRY_ALGORITHMS_HPP_
////////////////////////////////////////////////////////////////////////////////
namespace hom3 { namespace geometry {
////////////////////////////////////////////////////////////////////////////////

/// \brief Geometric algorithms
namespace algorithm {

template<class P1, class P2>
inline Num distance(const P1& p1, const P2& p2) { return (p1 - p2).norm(); }

template<int nd>
struct PointLT {
  template<class P1, class P2>
  inline constexpr bool operator()(const P1& a, const P2& b) const {
    if(a(0) > b(0)) { // \todo optimize?
      return false;
    } else if(math::approx(a(0),b(0))) {
      if(a(1) > b(1)) {
        return false;
      } else if(math::approx(a(1),b(1))) {
        if(nd == 2) {
          return false;
        } else if(nd == 3 && (a(2) > b(2) || math::approx(a(2),b(2)))) {
          return false;
        }
      }
    }
    return true;
  }
};

struct PointEQ {
  template<class P1, class P2>
  inline constexpr bool operator()(const P1& a, const P2& b) const {
    ASSERT(a.size() == b.size(), "Comparing vectors of different sizes!");
    const auto n = a.size();
    for(SInd d = 0; d < n; ++d) {
      if(!math::approx(a(d), b(d))) { return false; }
    }
    return true;
  }
};


template<SInd nd, class P1, class P2>
inline NumA<nd> cut_point_along_direction(P1 p1, const P2 p2, const Num lsv1,
                                          const Num lsv2, const SInd d) {
  const Num m = (lsv2 - lsv1) / (p2(d) - p1(d));
  const Num c = lsv1 - m * p1(d);
  const Num x = -c / m;
  p1(d) = x;
  return p1;
}

template<SInd nd, class Points, class Values>
inline Num linear_interpolation(NumA<nd> p, Points&& ps, Values&& vs) {
  const Num t = (p - ps(0)).norm() / (ps(1) - ps(0)).norm();
  return t * vs(0) + (1. - t) * vs(1);
}

template <SInd nd> struct surface_centroid {};

template <>
struct surface_centroid<2> {
  template<class Points>
  inline NumA<2> operator()(const Points& ps) {
    return ps(0) + 0.5 * (ps(0) + ps(1));
  }
};

template <>
struct surface_centroid<3> {
  template<class Points>
  inline NumA<3> operator()(const Points&) {
    TERMINATE("surface centroid for 3D not implemented yet");
  }
};

template<SInd nd>
inline Num area_triangle(NumA<nd> p0_, NumA<nd> p1_, NumA<nd> p2_) {

  const NumA<3> p0 = NumA<3>{p0_(0), p0_(1), (nd == 3) ? p0_(2) : 0.};
  const NumA<3> p1 = NumA<3>{p1_(0), p1_(1), (nd == 3) ? p1_(2) : 0.};
  const NumA<3> p2 = NumA<3>{p2_(0), p2_(1), (nd == 3) ? p2_(2) : 0.};

  const NumA<3> u = p1 - p0;
  const NumA<3> v = p2 - p0;
  return 0.5 * (u.cross(v)).norm();
}

template<SInd nd, class Points, class Values>
inline Num barycentric_interpolation(NumA<nd> p, Points&& ps, Values&& vs) {
    const Num A0 = area_triangle<nd>(p, ps(1), ps(2));
    const Num A1 = area_triangle<nd>(p, ps(0), ps(2));
    const Num A2 = area_triangle<nd>(p, ps(1), ps(0));
    return (A0 * vs(0) + A1 * vs(1) + A2 * vs(2)) / (A0 + A1 + A2);
}

template<SInd nd, class Points, class LT = PointLT<nd>, class EQ = PointEQ>
void remove_duplicate_points(Points& points, LT lt = PointLT<nd>{},
                             EQ eq = PointEQ{}) {
  boost::erase(points,
               boost::unique<boost::return_found_end>(
                   boost::sort(points, lt), eq));
}

}  // namespace algorithm

////////////////////////////////////////////////////////////////////////////////
}  // namespace geometry
}  // namespace hom3
////////////////////////////////////////////////////////////////////////////////
#endif
