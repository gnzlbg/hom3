#ifndef HOM3_MISC_RANGES_HPP_
#define HOM3_MISC_RANGES_HPP_
////////////////////////////////////////////////////////////////////////////////
/// \brief This file
////////////////////////////////////////////////////////////////////////////////
#include <boost/range.hpp>
#include <boost/range/counting_range.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range/algorithm_ext.hpp>
#include <boost/algorithm/cxx11/any_of.hpp>
#include <boost/algorithm/cxx11/all_of.hpp>
#include <boost/algorithm/cxx11/none_of.hpp>
#include <boost/algorithm/cxx11/copy_if.hpp>
#include <boost/algorithm/cxx11/copy_n.hpp>
#include <boost/algorithm/cxx11/iota.hpp>
#include <boost/algorithm/cxx11/one_of.hpp>
#include <boost/algorithm/cxx11/partition_point.hpp>
#include <boost/range/numeric.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/range/adaptors.hpp>
#include <boost/range/any_range.hpp>
#include <boost/range/join.hpp>
////////////////////////////////////////////////////////////////////////////////
namespace hom3 {
////////////////////////////////////////////////////////////////////////////////

/// \brief Range utilities and helper classes
namespace range {

/// \brief Implementation details of range
namespace detail {

/// \brief Range<T> is a counting range, and Range<T*> is an iterator range
template<class T> struct RangeImpl {
  using type = boost::iterator_range<boost::counting_iterator<T>>;
};
template<class T> struct RangeImpl<T*> {
  using type = boost::iterator_range<T*>;
};
}  // namespace detail

/// \brief Range
template<class T> using Range = typename detail::RangeImpl<T>::type;


template<class T> using AnyRange
= typename boost::any_range<T, boost::forward_traversal_tag, T, Int>;

template<class T, class R> const AnyRange<T> anyfy(const R& rng) {
  return {rng};
}

/// \brief Range filter
template<class T> using RangeFilter
= decltype(boost::adaptors::filtered(std::function<bool(T)>()));

namespace detail {

/// \brief Filtered range implementation (works on RangeImpl)
template<class T> struct FRangeImpl {
  using type = boost::filtered_range<std::function<bool(T)>, const Range<T> >;
};
template<class T> struct FRangeImpl<T*> {  // removes the pointer!
  using type = boost::filtered_range<std::function<bool(T)>, const Range<T*>>;
};

}  // namespace detail

/// \brief Filtered range
template<class T> using FRange  = typename detail::FRangeImpl<T>::type;

/// \brief Filtered^2 range
template<class T>
using F2Range = boost::filtered_range<std::function<bool(T)>, const FRange<T>>;

/// \brief Range transformer
template<class From, class To> using RangeTransformer
= decltype(boost::adaptors::transformed(std::function<To(From)>()));


// brief Range algorithms (i.e. all Boost.Range algorithms)
namespace algorithm {
// makes all Boost.Range algorithms available
using namespace boost::algorithm;
using boost::join;
using namespace boost::adaptors;
}  // namespace algorithm

}  // namespace range

// all range utilities are available in the hom3 namespace
using namespace range;

}  // namespace hom3
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
namespace boost { namespace range_detail {  // (Note: adding stuff to
                                            // range_detail is pretty ugly)
/// Negates range filters:
template <typename T> auto operator!(filter_holder<T> const& f)
-> decltype(adaptors::filtered(std::not1(f.val))) {
  return adaptors::filtered(std::not1(f.val));
}

}  // namespace range_detail

}  // namespace boost
////////////////////////////////////////////////////////////////////////////////
#endif
