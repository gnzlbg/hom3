#ifndef HOM3_SOLVERS_FV_CNS_QUANTITIES_HPP_
#define HOM3_SOLVERS_FV_CNS_QUANTITIES_HPP_
////////////////////////////////////////////////////////////////////////////////
/// \brief Defines the Navier-Stokes quantities.
////////////////////////////////////////////////////////////////////////////////
#include "globals.hpp"
#include "solver/fv/euler/quantities.hpp"
////////////////////////////////////////////////////////////////////////////////
namespace hom3 { namespace solver { namespace fv { namespace cns {
////////////////////////////////////////////////////////////////////////////////

/// \brief Physical quantities:
/// - Convert from dimensionless-to-dimensioned variables and vice-versa.
/// - Computes free-stream quantities.
struct Quantities : euler::Quantities {
  /// \brief Quantities constructor
  ///
  /// Requires properties:
  /// - everything required by euler::quantities
  /// - free-stream Reynoldsh number
  /// - stagnation Prandtl number
  /// - Sutherland's constant
  explicit Quantities(io::Properties properties) noexcept
    : euler::Quantities(properties)
    , ReInfty_(io::read<quantity::Dimensionless>(properties, "ReInfty"))
    , Pr0_(io::read<quantity::Dimensionless>(properties, "Pr0"))
    , S0_(io::read<quantity::Temperature>(properties, "Sutherland0"))
    , sOverTref_(S0_ / T0_)
    , sOverTrefP1_(sOverTref_ + 1.0)
    , Re0_(Re_infinity() * mu_infinity() / (rho_infinity() * u_infinity())) {
    std::cerr << "Free-stream variables: "
              << "ReInfty = " << ReInfty_ << " "
              << "muInfty = " << mu_infinity()
              << "\n";
    std::cerr << "Stagnation variables: "
              << "Pr0 = " << Pr0_ << " "
              << "S0 = " << S0_ << " "
              << "Re0 = " << Re0_ << " "
              << "mu0 = " << mu0()
              << "\n";
  }

  /// \name Input/Output Variables
  ///@{
  ///@}

  /// \name Dimensioned from dimensionless quantities
  ///@{

  /// \brief $\overline{\mu} = \mu \overline{\mu}_0$
  inline quantity::DynamicViscosity mu
  (const quantity::Dimensionless dynamicViscosity) const noexcept {
    return dynamicViscosity * mu0();
  }

  ///@}

  /// \name Dimensionless from dimensioned quantities
  ///@{

  /// \brief $\mu = \overline{\mu} / \overline{\mu}_0$
  inline quantity::Dimensionless mu
  (const quantity::DynamicViscosity dynamicViscosity) const noexcept {
    return dynamicViscosity / mu0();
  }

  ///@}

  /// \name Free-stream quantities
  ///@{

  /// Free-stream dynamic viscosity \mu_\infty = (T_\infty)^{3}{2} \frac{1 +
  /// S/T_{\mathrm{ref}}}{T_\infty + S/T_{\mathrm{ref}}}
  inline quantity::Dimensionless mu_infinity() const noexcept
  { return sutherlands_law(T_infinity()); }

  /// \brief Free stream Reynolds' number: \mathrm{Re}_\infty = \frac{
  /// \overline{\rho}_\infty \overline{u}_\infty \overline{L}_\mathrm{ref} }
  /// { \overline{\mu}_\infty }
  inline quantity::Dimensionless Re_infinity() const noexcept
  { return ReInfty_; }

  ///@}

  /// \name Stagnation quantities
  ///@{

  /// \brief Stagnation Prandtl number:
  /// \mathrm{Pr}_{0} = \frac{ \mu_0 c_{p,0} }{ \lambda_{0} }
  inline quantity::Dimensionless Pr0() const noexcept { return Pr0_; }

  /// \brief Stagnation Reynold's number: \mathrm{Re}_0 = \frac{
  /// \mathrm{Re}_\infty \mu_\infty }{ \rho_\infty u_\infty }
  inline quantity::Dimensionless Re0() const noexcept { return Re0_; }

  /// \brief Stagnation dynamic viscosity: \mu_0 = \frac{ \overline{a}_0
  /// \overline{L}_\mathrm{ref} \overline{\rho}_0 }{ \mathrm{Re}_0 }
  inline quantity::DynamicViscosity mu0() const noexcept
  { return a0() * Lref_ * rho0_ / Re0_; }

  ///@}


  /// \brief Computes viscosity according to Sutherland's law for the
  /// dimensionless temperature \param T
  ///
  /// Sutherland's law: \mu (T) = (T)^{3/2} * \frac{1 + S / T_{\mathrm{ref}}}{T
  /// + S / T_{\mathrm{ref}}}
  ///
  /// Here sOverTref_ = S / T_ref , and
  ///      sOverTrefP1_ = 1 + S / T_ref
  ///
  inline Num sutherlands_law(const quantity::Dimensionless T) const noexcept {
    return std::pow(T, 1.5) * sOverTrefP1_ / (T + sOverTref_);
  }

 private:
  /// Free stream Reynolds' number: \mathrm{Re}_\infty = \frac{
  /// \overline{\rho}_\infty \overline{u}_\infty \overline{L}_\mathrm{ref} }
  /// { \overline{\mu}_\infty }
  const quantity::Dimensionless ReInfty_;
  /// Stagnation Prandtl number:
  /// \mathrm{Pr}_{0} = \frac{ \mu_0 c_{p,0} }{ \lambda_{0} }
  const quantity::Dimensionless Pr0_;
  /// Sutherland's constant
  const quantity::Temperature S0_;
  /// Sutherland's constant over reference temperature: S / T_\mathrm{ref}
  const quantity::Dimensionless sOverTref_;
  /// Sutherland's constant over reference temperature plus one: S /
  /// T_\mathrm{ref} + 1
  const quantity::Dimensionless sOverTrefP1_;
  /// Stagnation Reynold's number: \mathrm{Re}_0 = \frac{ \mathrm{Re}_\infty
  /// \mu_\infty }{ \rho_\infty u_\infty }
  const quantity::Dimensionless Re0_;
};

/// \brief Computes conservative variables from primitive variables
template<SInd nd> static inline NumA<Indices<nd>::nvars> cv
(const NumA<Indices<nd>::nvars> pvars, const Num gamma)
{ return euler::cv<nd>(pvars, gamma); }

/// \brief Computes primitive variables from conservative variables
template<SInd nd> static inline NumA<Indices<nd>::nvars> pv
(const NumA<Indices<nd>::nvars> cvars, const Num gamma)
{ return euler::pv<nd>(cvars, gamma); }

////////////////////////////////////////////////////////////////////////////////
} // namespace cns
} // namespace fv
} // namespace solver
} // namespace hom3
////////////////////////////////////////////////////////////////////////////////
#endif
