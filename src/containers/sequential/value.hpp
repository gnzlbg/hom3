#ifndef HOM3_CONTAINER_SEQUENTIAL_VALUE_HPP_
#define HOM3_CONTAINER_SEQUENTIAL_VALUE_HPP_
////////////////////////////////////////////////////////////////////////////////
/// \file \brief Implementation of a container::sequential::value_type facade.
////////////////////////////////////////////////////////////////////////////////
/// Include:
#include <algorithm>
#include "traits.hpp"
////////////////////////////////////////////////////////////////////////////////
namespace hom3 { namespace container { namespace sequential {
////////////////////////////////////////////////////////////////////////////////

/// \brief Facade for a container::sequential::value_type
///
/// Value types are used to construct temporary container elements.
template<class Container> struct ValueFacade {
  /// \name Traits
  ///@{
  using value_type = typename traits<Container>::value_type;
  using reference  = typename traits<Container>::reference;
  ///@}

  /// \brief Underlying value type
  inline       value_type* r()       noexcept
  { return static_cast<value_type*>(this); }
  inline const value_type* r() const noexcept
  { return static_cast<const value_type*>(this); }

  friend inline void swap(value_type& lhs, value_type& rhs) noexcept
  { value_type::swap_values(lhs, rhs); }
};

////////////////////////////////////////////////////////////////////////////////
}  // namespace sequential
}  // namespace container
}  // namespace hom3
////////////////////////////////////////////////////////////////////////////////
#endif
