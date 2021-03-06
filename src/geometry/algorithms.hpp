#ifndef HOM3_GEOMETRY_ALGORITHMS_HPP_
#define HOM3_GEOMETRY_ALGORITHMS_HPP_
////////////////////////////////////////////////////////////////////////////////
namespace hom3 { namespace geometry {
////////////////////////////////////////////////////////////////////////////////

/// \brief Geometric algorithms
namespace algorithm {

template<class P1, class P2>
inline Num distance(const P1& p1, const P2& p2) { return (p1 - p2).norm(); }

}  // namespace algorithm

////////////////////////////////////////////////////////////////////////////////
}  // namespace geometry
}  // namespace hom3
////////////////////////////////////////////////////////////////////////////////
#endif
