#ifndef PROPERTY_INPUT_HPP_
#define PROPERTY_INPUT_HPP_
////////////////////////////////////////////////////////////////////////////////
/// Includes:
#include <utility>
#include <string>
#include <boost/any.hpp>
#include <unordered_map>
#include <type_traits>
////////////////////////////////////////////////////////////////////////////////
namespace hom3 {
////////////////////////////////////////////////////////////////////////////////

/// \brief Contains input/output functionality
namespace io {
////////////////////////////////////////////////////////////////////////////////

/// Property container
using Properties = std::unordered_map<std::string,boost::any>;
/// Property type
using Property = std::pair<std::string,boost::any>; ///< Property.

/// \brief Makes a property with \p name and \p value
Property make_property(std::string name, boost::any value) {
  return Property(name,value);
}

/// \brief Inserts a property with \p name and \p value into the
/// property container \p properties
///
/// \warning If a property with the same name already exists, this replaces its
/// value!
///
/// Note: std::remove_reference is used as an identity metafunction to prevent
/// template argument deduction.
template<class T> void insert_property (Properties& properties, std::string name,
typename std::remove_reference<T>::type value = typename std::remove_reference<T>::type()) {
  properties.insert(make_property(name,std::move(value)));
}

template<class T> void insert (Properties& properties, std::string name,
typename std::remove_reference<T>::type value = typename std::remove_reference<T>::type()) {
  properties.insert(make_property(name,std::move(value)));
}

/// \brief Reads property of type \p T with \p name from container \p
/// properties
///
/// \returns readed property of type \p T
template<class T>
T read(const Properties& properties, const std::string& name) {
  auto foundIt = properties.find(name);
  if(foundIt != std::end(properties) ){
    auto ret = foundIt->second;
    return boost::any_cast<T>(ret);
  } else {
    TERMINATE("property \"" + name + "\" not found!");
  }
}

/// \brief Reads property with \p name from property container \p properties
/// into the \p value.
template<class T>
void read(const Properties& properties, const std::string& name, T& value) {
  value = read<T>(properties, name);
}

////////////////////////////////////////////////////////////////////////////////
}} // hom3::io namespace
////////////////////////////////////////////////////////////////////////////////
#endif
