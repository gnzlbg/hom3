#ifndef HOM3_CONTAINER_SEQUENTIAL_ALGORITHMS_HPP_
#define HOM3_CONTAINER_SEQUENTIAL_ALGORITHMS_HPP_
////////////////////////////////////////////////////////////////////////////////
// Include:
#include <boost/range/algorithm.hpp>
#include <boost/algorithm/cxx11/all_of.hpp>
#include <boost/algorithm/cxx11/any_of.hpp>
#include <boost/algorithm/cxx11/copy_if.hpp>
#include <boost/algorithm/cxx11/copy_n.hpp>
#include <boost/algorithm/cxx11/iota.hpp>
#include <boost/algorithm/cxx11/none_of.hpp>
#include <boost/algorithm/cxx11/one_of.hpp>
#include <boost/algorithm/cxx11/partition_point.hpp>
////////////////////////////////////////////////////////////////////////////////
namespace container { namespace sequential {
////////////////////////////////////////////////////////////////////////////////

/// \brief Range algorithms that operate on cell containers
namespace algorithm {


/// Bring all boost::range::algorithms into the namespace
using namespace boost::algorithm;

/// \name Non-modifying operations
///@{

/// \brief Finds the first element in range [\p first, \p last) that
/// satisfies the predicate \p p
///
/// For containers: the predicate uses an index for accessing the container
/// _by reference_.
///
/// \algorithm non-modifying
/// \complexity O(n)
///@{
template <class T, class Predicate> T find_if(T first, const T last, Predicate&& p) {
  while(first != last && !p(first)) {
    ++first;
  }
  return first;
}

/// \brief Finds the first element in range \p c, that satisfies the predicate \p p
template<class Container, class Predicate>
auto find_if(Container&& c, Predicate&& p) {
  return find_if(c.first(),c.last(),std::forward<Predicate>(p));
}
///@}

///@}

/// \name Modifying operations
///@{

/// \brief Deletes elements in range \p c that satisfy the predicate \p p
///
/// \algorithm modifying
/// \complexity O(n) where n = #of elements in the range \p c
template<class Container, class Predicate>
Container& erase_remove_if(Container&& c, Predicate&& p) {
  c.pop_cell(std::end(std::forward<Container>(c))
             - boost::remove_if(std::forward<Container>(c),
                                std::forward<Predicate>(p)));
  // Maybe an indexed version is faster, e.g. something like
  // c.erase_remove_if(p);
  return std::forward<Container>(c);
}
///@}

////////////////////////////////////////////////////////////////////////////////
} } } // container::sequential::algorithm
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
#endif
