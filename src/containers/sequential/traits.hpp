#ifndef HOM3_CONTAINER_SEQUENTIAL_TRAITS_HPP_
#define HOM3_CONTAINER_SEQUENTIAL_TRAITS_HPP_
////////////////////////////////////////////////////////////////////////////////
/// \file \brief Trait class for sequential container.
////////////////////////////////////////////////////////////////////////////////
namespace container { namespace sequential {
////////////////////////////////////////////////////////////////////////////////

/// \brief Traits to distinguish between containers with fixed #of nodes per
/// element and those with variable #of nodes
template<class C> struct traits {};

/// \brief Tags to overload/specialize based on fixed/variable #of nodes per
/// container
namespace tag {
struct variable_nodes {};
struct fixed_nodes {};
}

////////////////////////////////////////////////////////////////////////////////
}} // container::sequential namespace
////////////////////////////////////////////////////////////////////////////////
#endif
