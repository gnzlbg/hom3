#ifndef HOM3_SOLVERS_FV_HEAT_QUANTITIES_HPP_
#define HOM3_SOLVERS_FV_HEAT_QUANTITIES_HPP_
////////////////////////////////////////////////////////////////////////////////
/// \brief Defines the Heat-equation's quantities.
////////////////////////////////////////////////////////////////////////////////
#include "globals.hpp"
////////////////////////////////////////////////////////////////////////////////
namespace hom3 { namespace solver { namespace fv { namespace heat {
////////////////////////////////////////////////////////////////////////////////

/// \brief Physical quantities:
/// - Convert from dimensionless-to-dimensioned variables and vice-versa.
struct Quantities {
  explicit Quantities(io::Properties properties) noexcept
    : alpha_(io::read<quantity::ThermalDiffusivity>(properties, "diffusivity"))
    , Tref_(io::read<quantity::Temperature>(properties, "T_ref"))
  {}

  /// \name Dimensioned from dimensionless quantities
  ///@{

  /// \brief $\overline{T} = T \overline{T}_0$
  inline quantity::Temperature T
  (const quantity::Dimensionless temperature) const noexcept
  { return temperature * Tref_; }

  /// \brief $\overline{t} = t \overline{L}_\mathrm{ref} / \overline{a}_0$
  inline quantity::Time t
  (const quantity::Dimensionless time) const noexcept
  { return time * unit::pow<2>(Lref_) / alpha_; }

  ///@}

  /// \name Dimensionless from dimensioned quantities
  ///@{

  /// \brief $T = \overline{T} / \overline{T}_0$
  inline quantity::Dimensionless T
  (const quantity::Temperature temperature) const noexcept
  { return temperature / Tref_; }

  /// \brief $t = \overline{t} \overline{a}_0  / \overline{L}_\mathrm{ref}$
  inline quantity::Dimensionless t
  (const quantity::Time time) const noexcept
  { return time * alpha_ / unit::pow<2>(Lref_); }

  ///@}

 private:
  /// Thermal diffusivity $\alpha = \lambda / (\rho c_p) \; [m^2 / s]$
  const quantity::ThermalDiffusivity alpha_;
  /// Reference temperature $T_\mathrm{ref} \; [K]$
  const quantity::Temperature Tref_;
  /// Reference length: $L_\mathrm{ref} \; [m]$
  const quantity::Length Lref_ {1. * unit::meter};
};

////////////////////////////////////////////////////////////////////////////////
} // namespace heat
} // namespace fv
} // namespace solver
} // namespace hom3
////////////////////////////////////////////////////////////////////////////////
#endif
