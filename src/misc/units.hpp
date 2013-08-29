#ifndef HOM3_MISC_UNITS_HPP_
#define HOM3_MISC_UNITS_HPP_
////////////////////////////////////////////////////////////////////////////////
/// \file \brief Includes unit/quantity types
////////////////////////////////////////////////////////////////////////////////
/// Includes:
#include "types.hpp"
#include <boost/units/systems/si.hpp>
#include <boost/units/cmath.hpp>
#include <boost/units/systems/temperature/celsius.hpp>
#include <boost/units/absolute.hpp>
#include <boost/units/detail/utility.hpp>
#include <boost/units/systems/si/io.hpp>
////////////////////////////////////////////////////////////////////////////////

namespace boost {

namespace units {

/// Derived dimension for energy density: L^-1, M, T^-2
using energy_density_dimension
  = derived_dimension<length_base_dimension, -1,
                      mass_base_dimension, 1,
                      time_base_dimension, -2>::type;

/// Derived dimension for specific energy: L^2, T^-2
using specific_energy_dimension
  = derived_dimension<length_base_dimension, 2,
                      time_base_dimension, -2>::type;

/// Derived dimension for the dynamic viscosity: M, L^-1, T^-1
using dynamic_viscosity_dimension
  = derived_dimension<mass_base_dimension, 1,
                      length_base_dimension, -1,
                      time_base_dimension, -1>::type;

/// Derived dimension for the thermal diffusivity: L^2, T^-1
using thermal_diffusivity_dimension
  = derived_dimension<length_base_dimension, 2,
                      time_base_dimension, -1>::type;

/// Derived dimension for the specific gas constant: L^2, T^-2, K^-1
using specific_gas_constant_dimension
  = derived_dimension<length_base_dimension, 2,
                      time_base_dimension, -2,
                      temperature_base_dimension, -1>::type;

namespace si {

using energy_density = unit<energy_density_dimension, si::system>;
BOOST_UNITS_STATIC_CONSTANT(joule_per_cubic_meter, energy_density);

using specific_energy = unit<specific_energy_dimension, si::system>;
BOOST_UNITS_STATIC_CONSTANT(joule_per_kilogram, specific_energy);

using dynamic_viscosity = unit<dynamic_viscosity_dimension, si::system>;
BOOST_UNITS_STATIC_CONSTANT(kg_per_meter_per_second, dynamic_viscosity);

using thermal_diffusivity = unit<thermal_diffusivity_dimension, si::system>;
BOOST_UNITS_STATIC_CONSTANT(meter_squared_per_second, thermal_diffusivity);

using specific_gas_constant = unit<specific_gas_constant_dimension, si::system>;
BOOST_UNITS_STATIC_CONSTANT(joule_per_kilogram_kelvin, specific_gas_constant);

} // namespace si

} // namespace units

} // namespace boos

namespace hom3 {

/// \brief Quantity Types
namespace quantity {

/// \brief Numeric quantity
template<class Dimensions>
using NumQ = boost::units::quantity<Dimensions, Num>;

using Dimensionless  = NumQ<boost::units::si::dimensionless>;    ///< [-]
using Length         = NumQ<boost::units::si::length>;           ///< [m]
using Time           = NumQ<boost::units::si::time>;             ///< [s]
using Mass           = NumQ<boost::units::si::mass>;             ///< [kg]
using Temperature    = NumQ<boost::units::si::temperature>;      ///< [K]
using Density        = NumQ<boost::units::si::mass_density>;     ///< [kg/m^3]
using Velocity       = NumQ<boost::units::si::velocity>;         ///< [m/s]
using Pressure       = NumQ<boost::units::si::pressure>;         ///< [kg/(ms^2)]
using EnergyDensity  = NumQ<boost::units::si::energy_density>;   ///< [J/m^3]
using SpecificEnergy = NumQ<boost::units::si::specific_energy>;  ///< [J/kg]

using DynamicViscosity
  = NumQ<boost::units::si::dynamic_viscosity>;  ///< [kg/(m s)]

using ThermalDiffusivity
  = NumQ<boost::units::si::thermal_diffusivity>;  ///< [m^2/s]

using SpecificGasConstant
  = NumQ<boost::units::si::specific_gas_constant>;  ///< [J/(kg K)]

} // namespace quantity

/// \brief Unit types
namespace unit {

using namespace boost::units;
using namespace boost::units::si;

} // namespace unit

////////////////////////////////////////////////////////////////////////////////
} // namespace hom3
////////////////////////////////////////////////////////////////////////////////
#endif
