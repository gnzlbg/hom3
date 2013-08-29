#ifndef HOM3_SOLVERS_FV_EULER_QUANTITIES_HPP_
#define HOM3_SOLVERS_FV_EULER_QUANTITIES_HPP_
////////////////////////////////////////////////////////////////////////////////
/// \brief Defines the Euler-equations' quantities.
////////////////////////////////////////////////////////////////////////////////
#include "globals.hpp"
////////////////////////////////////////////////////////////////////////////////
namespace hom3 { namespace solver { namespace fv { namespace euler {
////////////////////////////////////////////////////////////////////////////////

/// \brief Physical quantities:
/// - Convert from dimensionless-to-dimensioned variables and vice-versa.
/// - Computes free-stream quantities.
struct Quantities {

  /// \brief Quantities constructor
  ///
  /// Requires properties:
  /// - $\gamma$ (ratio of specific heats)
  /// - $R_\mathrm{specific} (specific gas constant)
  /// - $T_0$ (stagnation temperature)
  /// - $\rho_0$ (stagnation density)
  /// - $M_\infty$ (free-stream Mach number)
  explicit Quantities(io::Properties properties) noexcept
    : gamma_(io::read<quantity::Dimensionless>(properties, "gamma"))
    , Rspecific_
      (io::read<quantity::SpecificGasConstant>(properties, "Rspecific"))
    , T0_(io::read<quantity::Temperature>(properties, "T0"))
    , rho0_(io::read<quantity::Density>(properties, "rho0"))
    , MInfty_(io::read<quantity::Dimensionless>(properties, "Minfty"))
  {
    std::cerr << "Free-stream variables: "
              << "MInfty = " << MInfty_ << " "
              << "pInfty = " << p_infinity() << " "
              << "rhoInfty = " << rho_infinity() << " "
              << "TInfty = " << T_infinity() << " "
              << "uInfty = " << u_infinity()
              << "\n";
    std::cerr << "Stagnation variables: "
              << "a0 = " << a0() << " "
              << "rho0 = " << rho0_ << " "
              << "T0 = " << T0_
              << "\n";
    std::cerr << "Other variables: "
              << "gamma = " << gamma_ << " "
              << "Rspecific = " << Rspecific_ << " "
              << "Lref = " << Lref_ << " "
              << "\n";
  }

  /// \name Input/Output Variables
  ///@{

  /// \brief Stagnation speed of sound: $\overline{a}_0 = \sqrt{ \gamma
  /// \overline{R}_\mathrm{specific} \overline{T}_0 }$
  inline quantity::Velocity a0() const noexcept {
    return unit::sqrt(gamma() * Rspecific_ * T0_);
  }

  /// \brief Ratio of specific heats $\gamma = c_p / c_V$
  inline quantity::Dimensionless gamma() const noexcept { return gamma_; }
  /// \brief $\gamma - 1$
  inline quantity::Dimensionless gammaM1() const noexcept
  { return gamma() - 1.; }

  ///@}

  /// \name Dimensioned from dimensionless quantities
  ///@{

  /// \brief $\overline{p} = p \overline{\rho}_0 \overline{a}_0^2$
  inline quantity::Pressure p
  (const quantity::Dimensionless pressure) const noexcept
  { return pressure * rho0_ * unit::pow<2>(a0()); }

  /// \brief $\overline{\rho} = \rho \overline{\rho}_0$
  inline quantity::Density rho
  (const quantity::Dimensionless density) const noexcept
  { return density * rho0_; }

  /// \brief $\overline{u} = u \overline{a}_0$
  inline quantity::Velocity u
  (const quantity::Dimensionless velocity) const noexcept
  { return velocity * a0(); }

  /// \brief $\overline{t} = t \overline{L}_\mathrm{ref} / \overline{a}_0$
  inline quantity::Time t
  (const quantity::Dimensionless time) const noexcept
  { return time * Lref_ / a0() ; }

  /// \brief $\overline{T} = T \overline{T}_0$
  inline quantity::Temperature T
  (const quantity::Dimensionless temperature) const noexcept
  { return temperature * T0_; }
  ///@}

  /// \name Dimensionless from dimensioned quantities
  ///@{

  /// \brief $p = \overline{p} / ( \overline{\rho}_0 \overline{a}_0^2 )$
  inline quantity::Dimensionless p
  (const quantity::Pressure pressure) const noexcept
  { return pressure / (rho0_ * unit::pow<2>(a0())); }

  /// \brief $\rho = \overline{rho} / \overline{\rho}_0$
  inline quantity::Dimensionless rho
  (const quantity::Density density) const noexcept
  { return density / rho0_; }

  /// \brief $u = \overline{u} / \overline{a}_0$
  inline quantity::Dimensionless u
  (const quantity::Velocity velocity) const noexcept
  { return velocity / a0(); }

  /// \brief $t = \overline{t} \overline{a}_0  / \overline{L}_\mathrm{ref}$
  inline quantity::Dimensionless t
  (const quantity::Time time) const noexcept
  { return time * a0() / Lref_; }

  /// \brief $T = \overline{T} / \overline{T}_0$
  inline quantity::Dimensionless T
  (const quantity::Temperature temperature) const noexcept
  { return temperature / T0_; }
  ///@}

  /// \name Free-stream quantities
  ///@{

  /// \brief $M_\infty = \overline{u}_\infty / \overline{a}_\infty$
  inline quantity::Dimensionless M_infinity() const noexcept
  { return MInfty_; }

  /// \brief $T_\infty = (1 + (\gamma - 1) M_\infty^2 / 2)^{-1}$
  inline quantity::Dimensionless T_infinity() const noexcept
  { return unit::pow<-1>(1 + 0.5 * (gammaM1()) * unit::pow<2>(M_infinity())); }

  /// \brief $\rho_\infty = (T_\infty)^{1 / (\gamma - 1)}$
  inline quantity::Dimensionless rho_infinity() const noexcept
  { return unit::pow(T_infinity(), 1. / gammaM1()); }

  /// \brief $u_\infty = M_\infty \sqrt{T_\infty}$
  inline quantity::Dimensionless u_infinity() const noexcept
  { return M_infinity() * unit::sqrt(T_infinity()); }

  /// \brief $p_\infty = (T_\infty)^{ \gamma / (\gamma - 1) } / \gamma$
  inline quantity::Dimensionless p_infinity() const noexcept
  { return unit::pow(T_infinity(), gamma() / gammaM1()) / gamma(); }

  ///@}

  /// \name Primitive | Conservative variables conversions
  ///@{

  /// \brief Converts primitive variables \p pvs to conservative variables
  template<SInd nd> static inline NumA<Indices<nd>::nvars> cv
  (const NumA<Indices<nd>::nvars> pvs, const Num gamma) noexcept {
    using V =  Indices<nd>;
    NumA<V::nvars> cvs = NumA<V::nvars>::Zero();
    cvs(V::rho()) = pvs(V::rho());
    Num velMag2 = 0;
    for (SInd d = 0; d < nd; ++d) {
      cvs(V::rho_u(d)) = pvs(V::rho()) * pvs(V::u(d));
      velMag2 += std::pow(pvs(V::u(d)), 2);
    }
    cvs(V::rho_E()) = pvs(V::p()) / (gamma - 1.0)
                      + 0.5 * pvs(V::rho()) * velMag2;
    return cvs;
  }
  template<SInd nd> inline NumA<Indices<nd>::nvars> cv
  (const NumA<Indices<nd>::nvars> pvs) const noexcept
  { return cv<nd>(pvs, gamma()); }


  /// \brief Converts conservative variables \p pvs to primitive variables
  template<SInd nd> static inline NumA<Indices<nd>::nvars> pv
  (const NumA<Indices<nd>::nvars> cvs, const Num gamma) {
    using V =  Indices<nd>;
    NumA<V::nvars> pvs = NumA<V::nvars>::Zero();
    pvs(V::rho()) = cvs(V::rho());
    Num velMag2 = 0;
    for(SInd d = 0; d < nd; ++d) {
      pvs(V::u(d)) = cvs(V::rho_u(d)) / cvs(V::rho());
      velMag2 += std::pow(pvs(V::u(d)),2);
    }
    pvs(V::p()) = (gamma - 1.0) * (cvs(V::rho_E())
                                   - 0.5 * cvs(V::rho()) * velMag2);
    return pvs;
  }
  template<SInd nd> inline NumA<Indices<nd>::nvars> pv
  (const NumA<Indices<nd>::nvars> cvs) const noexcept
  { return pv<nd>(cvs, gamma()); }

  ///@}

 protected:
  /// Ratio of specific heats $\gamma = c_p / c_V \; [-]$
  const quantity::Dimensionless gamma_;
  /// Specific gas constant $R_\mathrm{specific} = R / m \; [J / (kg K)]$
  /// where $R \; [J / (mol K)]$ is the gas constant and $m \; [g / mol]$ is
  /// the molar mass of the gas mixture
  const quantity::SpecificGasConstant Rspecific_;
  /// Reference temperature $T_0 \; [K]$
  const quantity::Temperature T0_;
  /// Reference density $\rho_0 \; [kg / m^3]$
  const quantity::Density rho0_;
  /// Free-stream Mach number $M_\infty = u_\infty / a_\infty$
  const quantity::Dimensionless MInfty_;
  /// Reference length $L_\mathrm{ref} \; [m]$:
  const quantity::Length Lref_ { 1. * unit::meter };
};

/// \brief Computes conservative variables from primitive variables
template<SInd nd> static inline NumA<Indices<nd>::nvars> cv
(const NumA<Indices<nd>::nvars> pvs, const Num gamma) noexcept
{ return Quantities::template cv<nd>(pvs, gamma); }

/// \brief Computes primitive variables from conservative variables
template<SInd nd> static inline NumA<Indices<nd>::nvars> pv
(const NumA<Indices<nd>::nvars> cvs, const Num gamma)
{ return Quantities::template pv<nd>(cvs, gamma); }

////////////////////////////////////////////////////////////////////////////////
} // namespace euler
} // namespace fv
} // namespace solver
} // namespace hom3
////////////////////////////////////////////////////////////////////////////////
#endif
