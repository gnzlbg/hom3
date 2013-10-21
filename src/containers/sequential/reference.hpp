#ifndef HOM3_CONTAINER_SEQUENTIAL_REFERENCE_HPP_
#define HOM3_CONTAINER_SEQUENTIAL_REFERENCE_HPP_
////////////////////////////////////////////////////////////////////////////////
/// \file \brief Implementation of a container::sequential::reference proxy.
////////////////////////////////////////////////////////////////////////////////
/// Include:
#include <algorithm>
#include "traits.hpp"
////////////////////////////////////////////////////////////////////////////////
namespace hom3 { namespace container { namespace sequential {
////////////////////////////////////////////////////////////////////////////////

/// \brief Facade for a container::sequential::reference proxy
///
/// This proxy type provides the value_type interface for read/write acces to
/// both container elements and temporary value_types.
template<class Container> struct ReferenceFacade {
  /// \name Traits
  ///@{
  using value_type = typename traits<Container>::value_type;
  using reference  = typename traits<Container>::reference;
  using CIdx       = typename traits<Container>::cell_index_type;
  ///@}

  /// \brief Swap between two references
  friend inline void swap(reference   lhs, reference   rhs) noexcept
  { value_type::swap_values(lhs, rhs); }

  /// \brief Swap between a reference and a value
  friend inline void swap(reference   lhs, value_type& rhs) noexcept
  { value_type::swap_values(lhs, rhs); }
  friend inline void swap(value_type& lhs, reference   rhs) noexcept
  { swap(rhs, lhs); }

  /// \brief Returns the underlying reference type
  inline       reference* r()       noexcept
  { return static_cast<reference*>(this); }
  inline const reference* r() const noexcept
  { return static_cast<const reference*>(this); }

  inline bool is_reference_to_container() const noexcept
  { return c_ != nullptr; }

  /// \brief Assignment operator: from reference
  inline reference& operator=(reference rhs) noexcept {
    if (c() && c() == rhs.c()) {
      // both are references to the same container
      r()->c()->copy_cell(rhs.index(), index());
    } else {
      // one of the following:
      // - references to different containers
      // - this refers to a container and rhs to a value
      // - this refers to a value and rhs to a container
      // - both are references to values
      value_type::copy_values(rhs, *r());
    }
    return *r();
  }

  /// \brief Assignment operator: from value
  inline reference& operator=(value_type rhs) noexcept {
    value_type::copy_values(rhs, *r());
    return *r();
  }

  /// \brief Implicit conversion its value type
  operator value_type() const noexcept {
    value_type ret;
    value_type::copy_values(*r(), ret);
    return ret;
  }
  operator value_type() noexcept {
    value_type ret;
    value_type::copy_values(*r(), ret);
    return ret;
  }

  /// \brief Returns a pointer to the container being referenced
  /// \warning The pointer can be a nullptr !
  inline Container*  c() const noexcept { return c_; }
  inline Container*& c()       noexcept { return c_; }

  /// \brief Returns the index of the element within a container, if the
  /// reference referes to an element within a container
  inline CIdx  index() const noexcept {
    ASSERT(is_valid(index_), "ERROR!");
    return index_;
  }
  inline CIdx& index()       noexcept {
    ASSERT(is_valid(index_), "ERROR!");
    return index_;
  }

  ReferenceFacade() noexcept : c_(nullptr), index_(invalid<CIdx>()) {}
  ReferenceFacade(Container* c, CIdx i) noexcept : c_(c), index_(i) {}

 private:
  Container* c_;  ///< Pointer to the underlying container
  CIdx   index_;  ///< Index of the of the element within the container
};

////////////////////////////////////////////////////////////////////////////////
}  // namespace sequential
}  // namespace container
}  // namespace hom3
////////////////////////////////////////////////////////////////////////////////
#endif
