#ifndef HOM3_CONTAINER_SEQUENTIAL_ITERATOR_HPP_
#define HOM3_CONTAINER_SEQUENTIAL_ITERATOR_HPP_
////////////////////////////////////////////////////////////////////////////////
/// \file \brief This file contains the implementation of the container
/// iterator as well as the value and reference proxy types facades.
////////////////////////////////////////////////////////////////////////////////
/// Include:
#include "traits.hpp"
////////////////////////////////////////////////////////////////////////////////
namespace container { namespace sequential {
////////////////////////////////////////////////////////////////////////////////

#define assert_in_cell_range(cIdx)                                      \
  ASSERT(cIdx() <= c()->last(),                                         \
         "Iterator's cell index " << cIdx()                             \
         << " is out of bounds [" << c()->first() << ","                \
         << c()->size() <<  "].")

#define assert_same_containers(lhs,rhs)                         \
  ASSERT(same_container(lhs,rhs), "Different containers")

/// \brief Boilerplate for implementing a value proxy type
///
/// Value proxy types provide a value-like interface for the container
/// elements which is used by e.g. algorithms to construct temporary elements.
///
template<class Container> struct ValueFacade {
  /// \name Traits
  ///@{
  using value_type = typename traits<Container>::value_type;
  using reference  = typename traits<Container>::reference;
  ///@}

  /// \brief Returns the underlying value type
  value_type* r(){return static_cast<value_type*>(this);}
  const value_type* r() const {return static_cast<const value_type*>(this);}

  /// \brief Implicit conversion to its reference type
  operator reference() const {
    reference ref;
    value_type::copy_values(ref,const_cast<value_type&>(*r()));
    return ref;
  }
  operator reference() {
    reference ref;
    value_type::copy_values(ref,const_cast<value_type&>(*r()));
    return ref;
  }

  friend inline void swap(value_type& lhs, value_type& rhs) {
    value_type::swap_values(lhs,rhs);
  }
};

/// \brief Boilerplate for implementing a reference proxy type
///
/// Reference proxy types provide access to the container elements which is
/// used when iterators are dereferenced.
///
template<class Container> struct ReferenceFacade {
  /// \name Traits
  ///@{
  using value_type = typename traits<Container>::value_type;
  using reference  = typename traits<Container>::reference;
  using CIdx       = typename traits<Container>::cell_index_type;
  ///@}

  friend inline void swap(reference   lhs, reference   rhs) { value_type::swap_values(lhs,rhs); }
  friend inline void swap(reference   lhs, value_type& rhs) { value_type::swap_values(lhs,rhs); }
  friend inline void swap(value_type& lhs, reference   rhs) { swap(rhs,lhs); }

  /// \brief Returns the underlying reference type
  reference* r(){return static_cast<reference*>(this);}
  const reference* r() const {return static_cast<const reference*>(this);}

  /// \brief Assignment operator
  reference& operator=(reference rhs) {
    ASSERT(c_ == rhs.c(), "ERROR!");
    r()->c_->copy_cell(rhs.index(),index());
    return *r();
  }
  reference& operator=(value_type rhs) {
    value_type::copy_values(*r(),rhs);
    return *r();
  }

  /// \brief Implicit conversion its value type
  operator value_type() const {
    value_type ret;
    value_type::copy_values(ret,*r());
    return ret;
  }
  operator value_type() {
    value_type ret;
    value_type::copy_values(ret,*r());
    return ret;
  }

  /// \brief Returns a pointer to the container being referenced
  /// \warning The pointer can be a nullptr !
  inline Container*  c() const { ASSERT( c_ != nullptr, "ERROR!"); return c_; }
  inline Container*& c()       { ASSERT( c_ != nullptr, "ERROR!"); return c_; }

  /// \brief Returns the index of the element within a container, if the reference
  /// referes to an element within a container
  inline CIdx  index() const { ASSERT( is_valid(index_), "ERROR!"); return index_; }
  inline CIdx& index()       { ASSERT( is_valid(index_), "ERROR!"); return index_; }


  ReferenceFacade() : c_(nullptr), index_(invalid<CIdx>())  {}
  ReferenceFacade(Container* c, CIdx i) : c_(c), index_(i) {}
 private:
  Container* c_; ///< Pointer to the underlying container
  CIdx   index_; ///< Index of the of the element
};

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

  Iterator() : c_(), index_()  {}
  Iterator(Container* c, CIdx index) : c_(c), index_(index) {}

  /// \brief Returns a pointer to the iterators container
  inline Container* c()             { return c_;     }
  inline const Container* c() const { return c_;     }

  /// \brief Returns the index of the element being pointer
  /// by the iterator within the container
  inline CIdx index()            { return index_; }
  inline CIdx index()      const { return index_; }

  /// \name Comparison operators (==, !=, <, >, <=, >=)
  ///@{
  static inline bool same_container(const Iterator& lhs, const Iterator& rhs) {
    return lhs.c() == rhs.c();
  }
  friend inline bool operator==(const Iterator& lhs, const Iterator& rhs) {
    assert_same_containers(lhs,rhs);
    return lhs.index() == rhs.index();
  }
  friend inline bool operator<=(const Iterator& lhs, const Iterator& rhs) {
    assert_same_containers(lhs,rhs);
    return lhs.index() <= rhs.index();
  }
  friend inline bool operator>=(const Iterator& lhs, const Iterator& rhs) {
    assert_same_containers(lhs,rhs);
    return lhs.index() >= rhs.index();
  }
  friend inline bool operator!=(const Iterator& lhs, const Iterator& rhs) {
    return !(lhs == rhs);
  }
  friend inline bool operator<(const Iterator& lhs, const Iterator& rhs) {
    return !(lhs >= rhs);
  }
  friend inline bool operator>(const Iterator& lhs, const Iterator& rhs) {
    return !(lhs <= rhs);
  }
  ///@}

  /// \name Traversal operators (++, +=, +, --, -=, -)
  ///@{
  inline Iterator& operator++() {
    ++index_;
    assert_in_cell_range(index);
    return *this;
  }
  inline Iterator operator++(int) { return Iterator(c_, index_++); }
  inline Iterator& operator+=(const difference_type o) {
    primitive_cast(index_) += o;
    assert_in_cell_range(index);
    return *this;
  }
  friend inline Iterator operator+(Iterator it, difference_type o) {
    return it += o;
  }
  friend difference_type operator+(const Iterator& lhs, const Iterator& rhs) {
    assert_same_containers(lhs,rhs);
    return lhs.index_ + rhs.index_;
  }
  inline Iterator& operator--() {
    --index_;
    assert_in_cell_range(index);
    return *this;
  }
  inline Iterator operator--(int) { return Iterator(c_, index_--); }
  inline Iterator& operator-=(difference_type o) {
    primitive_cast(index_) -= o;
    assert_in_cell_range(index);
    return *this;
  }
  friend inline Iterator operator-(Iterator it, difference_type o) { return it -= o; }
  friend inline difference_type operator-(const Iterator& lhs, const Iterator& rhs) {
    assert_same_containers(lhs,rhs);
    return static_cast<difference_type>(primitive_cast(lhs.index_))
        -
        static_cast<difference_type>(primitive_cast(rhs.index_));
  }
  ///@}

  /// \name Access operators
  ///@{
  inline reference operator*()  const { return {c_, index_}; }
  inline reference operator->() const { return {c_, index_}; }
  inline reference operator[](const difference_type o) const { return {c_, index_ + o}; }
  ///@}

  friend inline void swap(Iterator& lhs, Iterator& rhs) {
    using std::swap;
    swap(lhs.c_, rhs.c_);
    swap(lhs.index_, rhs.index_);
  }

 private:
  Container* c_; ///< Pointer to the container
  CIdx index_;    ///< Index of the element
};

////////////////////////////////////////////////////////////////////////////////
}} // container::sequential namespace

////////////////////////////////////////////////////////////////////////////////
#undef assert_in_cell_range
#undef assert_same_containers
////////////////////////////////////////////////////////////////////////////////
#endif
