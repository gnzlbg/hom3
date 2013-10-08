#ifndef HOM3_CONTAINER_SEQUENTIAL_TRAITS_HPP_
#define HOM3_CONTAINER_SEQUENTIAL_TRAITS_HPP_
////////////////////////////////////////////////////////////////////////////////
/// \file \brief Trait class for sequential container.
////////////////////////////////////////////////////////////////////////////////
namespace hom3 { namespace container { namespace sequential {
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

/// \brief Type-Trait Introspection
namespace tti {

/// \brief Has \p T a node index type?
struct has_node_index_type {
  template<class T>
  static auto test() -> decltype(T::node_index_type(), std::true_type());
  template<class>
  static auto test(...) -> std::false_type;
};

template<class T, class D> constexpr auto get_node_idx_type(std::true_type)
{ return T::node_index_type(); }
template<class T, class D> constexpr auto get_node_idx_type(std::false_type)
{ return D(); }

/// \brief Returns the node index type from \p Trait or \p Default if
/// \p Trait has no node index type.
template<class Trait, class Default = SInd> constexpr auto node_index_type() {
  return get_node_idx_type<Trait, Default>(has_node_index_type::test<Trait>());
}

}  // namespace tti

////////////////////////////////////////////////////////////////////////////////
}  // namespace sequential
}  // namespace container
}  // namespace hom3
////////////////////////////////////////////////////////////////////////////////
#endif
