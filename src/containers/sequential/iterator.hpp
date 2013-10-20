#ifndef HOM3_CONTAINER_SEQUENTIAL_ITERATOR_HPP_
#define HOM3_CONTAINER_SEQUENTIAL_ITERATOR_HPP_
////////////////////////////////////////////////////////////////////////////////
/// \file \brief This file implements the container::sequential iterator type.
////////////////////////////////////////////////////////////////////////////////
/// Include:
#include <algorithm>
#include "containers/sequential/traits.hpp"
#include "containers/sequential/value.hpp"
#include "containers/sequential/reference.hpp"
////////////////////////////////////////////////////////////////////////////////
/// File macros:
#define assert_in_cell_range(cIdx)                                      \
  ASSERT(cIdx() <= c()->last(),                                         \
         "Iterator's cell index " << cIdx()                             \
         << " is out of bounds [" << c()->first() << ","                \
         << c()->size() <<  "].")

#define assert_same_containers(lhs, rhs)                        \
  ASSERT(same_container(lhs, rhs), "Different containers")
////////////////////////////////////////////////////////////////////////////////
namespace hom3 { namespace container { namespace sequential {
////////////////////////////////////////////////////////////////////////////////

/// \brief Random Access iterator implementation
///
/// This class provides a random access iterator. Each iterator contains:
///  - a pointer to the container
///  - an index to the element being accessed
///
/// The existence of elements is emulated with the value/reference types
/// \see ValueFacade, ReferenceFacade
template<class Container> struct Iterator {
  /// \name Traits
  ///@{
  using difference_type   = Int;
  using size_type         = Ind;
  using value_type        = typename traits<Container>::value_type;
  using reference         = typename traits<Container>::reference;
  using pointer           = typename traits<Container>::reference;
  using iterator_category = std::random_access_iterator_tag;
  using CIdx              = typename traits<Container>::cell_index_type;
  ///@}

  Iterator() noexcept : c_(), index_()  {}
  Iterator(Container* c, CIdx index) noexcept : c_(c), index_(index) {}

  /// \brief Returns a pointer to the iterators container
  inline       Container* c()       noexcept { return c_; }
  inline const Container* c() const noexcept { return c_; }

  /// \brief Returns the index of the element being pointer
  /// by the iterator within the container
  inline CIdx index()       noexcept { return index_; }
  inline CIdx index() const noexcept { return index_; }

  /// \name Comparison operators (==, !=, <, >, <=, >=)
  ///@{
  static inline
  bool same_container(const Iterator& lhs, const Iterator& rhs) noexcept
  { return lhs.c() == rhs.c(); }
  friend inline
  bool operator==(const Iterator& lhs, const Iterator& rhs) noexcept
  { assert_same_containers(lhs, rhs); return lhs.index() == rhs.index(); }
  friend inline
  bool operator<=(const Iterator& lhs, const Iterator& rhs) noexcept
  { assert_same_containers(lhs, rhs); return lhs.index() <= rhs.index(); }
  friend inline
  bool operator>=(const Iterator& lhs, const Iterator& rhs) noexcept
  { assert_same_containers(lhs, rhs); return lhs.index() >= rhs.index(); }
  friend inline
  bool operator!=(const Iterator& lhs, const Iterator& rhs) noexcept
  { return !(lhs == rhs); }
  friend inline
  bool operator<(const Iterator& lhs, const Iterator& rhs) noexcept
  { return !(lhs >= rhs); }
  friend inline
  bool operator>(const Iterator& lhs, const Iterator& rhs) noexcept
  { return !(lhs <= rhs); }
  ///@}

  /// \name Traversal operators (++, +=, +, --, -=, -)
  ///@{
  inline Iterator& operator++() noexcept {
    ++index_;
    assert_in_cell_range(index);
    return *this;
  }
  inline Iterator operator++(int) noexcept { return Iterator(c_, index_++); }
  inline Iterator& operator+=(const difference_type o) noexcept {
    primitive_cast(index_) += o;
    assert_in_cell_range(index);
    return *this;
  }
  friend inline Iterator operator+(Iterator it, difference_type o) noexcept
  { return it += o; }
  friend inline
  difference_type operator+(const Iterator& lhs, const Iterator& rhs) noexcept {
    assert_same_containers(lhs, rhs);
    return lhs.index_ + rhs.index_;
  }
  inline Iterator& operator--() noexcept {
    --index_;
    assert_in_cell_range(index);
    return *this;
  }
  inline Iterator operator--(int) noexcept { return Iterator(c_, index_--); }
  inline Iterator& operator-=(difference_type o) noexcept {
    primitive_cast(index_) -= o;
    assert_in_cell_range(index);
    return *this;
  }
  friend inline Iterator operator-(Iterator it, difference_type o) noexcept
  { return it -= o; }
  friend inline
  difference_type operator-(const Iterator& lhs, const Iterator& rhs) noexcept {
    assert_same_containers(lhs, rhs);
    return static_cast<difference_type>(primitive_cast(lhs.index_))
        - static_cast<difference_type>(primitive_cast(rhs.index_));
  }
  ///@}

  /// \name Access operators
  ///@{
  inline reference operator*()  const noexcept { return {c_, index_}; }
  inline reference operator->() const noexcept { return {c_, index_}; }
  inline reference operator[](const difference_type o) const noexcept
  { return {c_, index_ + o}; }
  ///@}

  friend inline void swap(Iterator& lhs, Iterator& rhs) noexcept {
    using std::swap;
    swap(lhs.c_, rhs.c_);
    swap(lhs.index_, rhs.index_);
  }

 private:
  Container* c_;  ///< Pointer to the container
  CIdx index_;    ///< Index of the element
};

////////////////////////////////////////////////////////////////////////////////
}  // namespace sequential
}  // namespace container
}  // namespace hom3
////////////////////////////////////////////////////////////////////////////////
#undef assert_in_cell_range
#undef assert_same_containers
////////////////////////////////////////////////////////////////////////////////
#endif
