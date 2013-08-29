#ifndef HOM3_SOLVERS_FV_CNS_PHYSICS_HPP_
#define HOM3_SOLVERS_FV_CNS_PHYSICS_HPP_
////////////////////////////////////////////////////////////////////////////////
/// \file \brief Implements the physics class of the Navier-Stokes equations
/// for ideal compressible flow.
////////////////////////////////////////////////////////////////////////////////
#include <limits>
#include <string>
#include <algorithm>
#include "solver/fv/solver.hpp"
#include "solver/fv/euler/physics.hpp"
#include "indices.hpp"
#include "tags.hpp"
#include "quantities.hpp"
/// Options:
#define ENABLE_DBG_ 0
#include "misc/dbg.hpp"
////////////////////////////////////////////////////////////////////////////////
namespace hom3 { namespace solver { namespace fv {
////////////////////////////////////////////////////////////////////////////////

/// \brief Navier-Stokes equations for ideal compressible flow
namespace cns {

/// \brief Physics class of the Navier-Stokes equations for ideal compressible
/// flow.
template<SInd nd_, class NumFlux, class Solver>
struct Physics : euler::Physics<nd_, NumFlux, Solver> {
  /// \name Solver interface
  ///@{
  using physics_type = type_tag;
  static const constexpr SInd nd = nd_;
  static std::string physics_name() noexcept {
    return "CNS" + std::to_string(nd) + "D";
  }
  ///@}

  /// \name Inherit behavior from Euler-eqts:
  ///@{
  using Euler = euler::Physics<nd_, NumFlux, Solver>;
  using Euler::b_; using Euler::gamma; using Euler::gammaM1;
  template<class _> inline Num& rho(const CellIdx cIdx)       noexcept
  { return Euler::template rho<_>(cIdx); }
  template<class _> inline Num  rho(const CellIdx cIdx) const noexcept
  { return Euler::template rho<_>(cIdx); }
  template<class _> inline Num  rho_u(const CellIdx cIdx, SInd d) const noexcept
  { return Euler::template rho_u<_>(cIdx, d); }
  template<class _> inline Num u
  (const CellIdx cIdx, const SInd d) const noexcept
  { return Euler::template u<_>(cIdx, d); }
  template<class _> inline NumA<nd> u(const CellIdx cIdx) const noexcept
  { return Euler::template u<_>(cIdx); }
  template<class _> inline Num p(const CellIdx cIdx) const noexcept
  { return Euler::template p<_>(cIdx); }
  ///@}

  /// Constructor: requires property gamma
  explicit Physics(io::Properties properties) noexcept
      : euler::Physics<nd_, NumFlux, Solver>(properties)
      , quantities(properties)
      , cfl_viscous_(io::read<Num>(properties, "CFL_viscous"))
  {}

  /// \name Variable Access
  ///@{

  /// Alias for data Indices
  using V = Indices<nd>;

  /// #of variables: nd + 2 = nd (u_vector) + 1 (rho) + 1 (E)
  static constexpr SInd nvars = V::nvars;

  /// \brief Derivative of the pressure in direction \p slopeDir at \p cIdx
  template<class _>
  Num slope_p(const CellIdx cIdx, const SInd slopeDir) const noexcept {
    return b_()->slope(cIdx, [&](const CellIdx i) { return p<_>(i); },
                       slopeDir);
  }

  /// \brief Derivative of the density in direction \p slopeDir at \p cIdx
  template<class _>
  Num slope_rho(const CellIdx cIdx, const SInd slopeDir) const noexcept {
    return b_()->slope(cIdx, [&](const CellIdx i) { return rho<_>(i); },
                       slopeDir);
  }

  /// \brief Derivative of the \p dir velocity component in direction \p
  /// slopeDir at \p cIdx
  template<class _>
  Num slope_u(const CellIdx cIdx, const SInd dir,
              const SInd slopeDir) const noexcept {
    return (b_()->template slope<_>(cIdx, V::rho_u(dir), slopeDir) * rho<_>(cIdx)
            - b_()->template slope<_>(cIdx, V::rho(), slopeDir) * rho_u<_>(cIdx, dir))
        / std::pow(rho<_>(cIdx), 2.);
  }

  /// \brief Computes the heat flux
  ///
  /// Heat flux:
  /// q = - \frac{\mu \gamma}{(\gamma - 1) \mathrm{Pr}_0 \rho}
  /// \left ( \nabla p - \frac{p}{\rho} \nabla \rho \right )
  inline Num heat_flux(const Num rho, const Num gradRho, const Num mu,
                       const Num p, const Num gradP) const noexcept {
    return - mu * gamma() / (gammaM1() * quantities.Pr0() * rho)
        * (gradP - p / rho * gradRho);
  }

  /// \brief Computes the shear stress tensor
  /// \tau = \mu (- (\nabla \mathbf{u} + (\nabla \mathbf{u})^{T})
  ///             + \frac{2}{3} (\nabla \cdot \mathbf{u}) \mathbf{I})
  template<class S>
  NumAM<nd, nd> shear_stress(const Num mu, S&& slope_u) const noexcept {
    NumAM<nd, nd> nabla_u; // use binaryExpr instead!
    for(SInd row_i = 0; row_i < nd; ++row_i) {
      for(SInd col_j = 0; col_j < nd; ++col_j) {
        nabla_u(row_i, col_j) = slope_u(col_j, row_i);
      }
    }
    DBGV((mu)(nabla_u));
    const Num div_u = nabla_u.diagonal().sum() * 2. / 3.;

    return - mu * (nabla_u + nabla_u.transpose()
                   - div_u * NumAM<nd, nd>::Identity());
  }

  template<class _>
  inline NumA<1> vorticity(const CellIdx cIdx, grid::dim<2>) const noexcept {
    NumA<1> vorticity;
    vorticity(0) = slope_u<_>(cIdx, 1, 0) - slope_u<_>(cIdx, 0, 1);
    return vorticity;
  }

  template<class _>
  inline NumA<nd> vorticity(const CellIdx cIdx, grid::dim<3>) const noexcept {
    NumA<nd> vorticity = NumA<nd>::Zero();
    vorticity(0) = slope_u<_>(cIdx, 2, 1) - slope_u<_>(cIdx, 1, 2);
    vorticity(1) = slope_u<_>(cIdx, 0, 2) - slope_u<_>(cIdx, 2, 0);
    vorticity(2) = slope_u<_>(cIdx, 1, 0) - slope_u<_>(cIdx, 0, 1);
  }

  template<class _>
  inline auto vorticity(const CellIdx cIdx) const
  RETURNS(vorticity<_>(cIdx, grid::dim<nd>()));

  static inline constexpr SInd no_vorticity_components() noexcept
  { return nd == 2 ? 1 : 3; }

  ///@}

  /// \name Numerical functions
  ///@{

  /// \brief Computes the numerical flux at the interface between \p lIdx and
  /// \p rIdx
  template<class _>
  inline NumA<nvars> compute_num_flux
  (const CellIdx lIdx, const CellIdx rIdx, const SInd d, const Num dx,
   const Num dt) const noexcept {
    return Euler::template
        compute_num_flux_<_>(lIdx, rIdx, d, dx, dt, NumFlux())
        + compute_num_flux_viscous_<_>(lIdx, rIdx, d, dx, dt, NumFlux());
  }

  /// \brief computes dt at cell \p lIdx
  ///
  /// $\min_{dt_\mathrm{inviscid}, dt_\mathrm{viscous}}$ where $dt_viscous =
  /// \frac{h^2}{ \mu / \rho }$
  template<class _> inline Num compute_dt(const CellIdx cIdx) const noexcept {
    const auto dt_inviscid = Euler::template compute_dt<_>(cIdx);

    const Num h = b_()->cells().length(cIdx);

    /// Compute viscous "wave speed": S_viscous = \mu / (\rho \delta x)
    const Num S_coeff = mu<_>(cIdx) /  (rho<_>(cIdx) * h);
    /// Compute viscous coefficient for shear-stress and heat conduction terms:
    /// - shear stress: S_coeff / Re_0
    /// - heat conduction: S_coeff / (Re_0 Pr_0 (\gamma - 1))
    const Num S_viscous = std::max
      (S_coeff /  quantities.Re0(),
       S_coeff / (quantities.Re0() * quantities.Pr0() * quantities.gammaM1()));

    const Num dt_viscous = cfl_viscous_ * h  / (2 * S_viscous);

    const Num dt = std::min(dt_inviscid, dt_viscous);

    DBGV((cIdx)(dt_inviscid)(dt_viscous)(dt)(S_viscous)(cfl_viscous_));
    return dt;
  }

  ///@}

  /// \todo add shear stress
  /// \todo add stress tensor
  template<class Output> void physics_output(Output&& out) const noexcept {
    Euler::physics_output(out);
    out << io::stream("mu", 1, [&](const Ind cIdx, const SInd) {
      return mu<lhs_tag>(CellIdx{cIdx});
    });
    out << io::stream("vorticity", no_vorticity_components(),
                      [&](const Ind cIdx, const SInd d)
                      { return vorticity<lhs_tag>(CellIdx{cIdx})(d);});
  }

  /// Physical quantities
  const Quantities quantities;

 private:
  const Num cfl_viscous_;

  static inline Num temperature
  (const Num gamma, const Num p, const Num rho) noexcept
  { return gamma * p / rho; }
  inline Num mu(const Num T) const noexcept
  { return quantities.sutherlands_law(T); }
  template<class _> inline Num mu(const CellIdx cIdx) const noexcept
  { return quantities.sutherlands_law(this->template T<_>(cIdx)); }

  /// Numerical fluxes implementation
  ///@{

  /// \todo this functions are pretty ugly: clean them up
  template<class _>
  inline Num rho_srfc(const CellIdx lIdx, const CellIdx rIdx) const noexcept
  { return b_()->Q_srfc(_(), lIdx, rIdx, V::rho()); }

  template<class _>
  inline Num rho_u_srfc(const CellIdx lIdx, const CellIdx rIdx,
                        const SInd d) const noexcept
  { return b_()->Q_srfc(_(), lIdx, rIdx, V::rho_u(d)); }

  template<class _>
  inline Num rho_E_srfc(const CellIdx lIdx, const CellIdx rIdx) const noexcept
  { return b_()->Q_srfc(_(), lIdx, rIdx, V::rho_E()); }

  template<class _>
  inline Num u_srfc(const CellIdx lIdx, const CellIdx rIdx,
                    const SInd d) const noexcept {
    return rho_u_srfc<_>(lIdx, rIdx, d) / rho_srfc<_>(lIdx, rIdx);
  }

  template<class _>
  inline Num umag2_srfc(const CellIdx lIdx, const CellIdx rIdx) const noexcept {
    Num umag2 = 0;
    for(SInd d = 0; d < nd; ++d) {
      umag2 += std::pow(u_srfc<_>(lIdx, rIdx, d), 2.);
    }
    return umag2;
  }

  template<class _>
  inline Num p_srfc(const CellIdx lIdx, const CellIdx rIdx) const noexcept {
    return gammaM1() * (rho_E_srfc<_>(lIdx, rIdx)
                        - 0.5 * rho_srfc<_>(lIdx, rIdx)
                        * umag2_srfc<_>(lIdx, rIdx));
  }

  template<class _>
  inline Num grad_rho_srfc(const CellIdx lIdx, const CellIdx rIdx,
                           const SInd d) const noexcept
  { return b_()->surface_slope(_(), lIdx, rIdx, V::rho(), d); }

  template<class _>
  inline Num grad_rho_u_srfc(const CellIdx lIdx, const CellIdx rIdx,
                             const SInd c, const SInd d) const noexcept
  { return b_()->surface_slope(_(), lIdx, rIdx, V::rho_u(c), d); }

  template<class _>
  inline Num grad_rho_E_srfc(const CellIdx lIdx, const CellIdx rIdx,
                             const SInd d) const noexcept
  { return b_()->surface_slope(_(), lIdx, rIdx, V::rho_E(), d); }


  template<class _>
  inline Num grad_u_srfc(const CellIdx lIdx, const CellIdx rIdx,
                         const SInd c, const SInd d) const noexcept {
    const auto rho_s = rho_srfc<_>(lIdx, rIdx);
    return (grad_rho_u_srfc<_>(lIdx, rIdx, c, d) * rho_s
            - grad_rho_srfc<_>(lIdx, rIdx, d) * rho_u_srfc<_>(lIdx, rIdx, c))
        / std::pow(rho_s, 2.);
  }

  template<class _>
  inline Num grad_p_srfc(const CellIdx lIdx, const CellIdx rIdx,
                         const SInd dir) const noexcept {
    const Num rho_srfc_ = rho_srfc<_>(lIdx, rIdx);
    Num grad_rho_umag2 = grad_rho_srfc<_>(lIdx, rIdx, dir)
                         * umag2_srfc<_>(lIdx, rIdx);
    for(SInd d = 0; d < nd; ++d) {
      grad_rho_umag2 += rho_srfc_ * 2. * u_srfc<_>(lIdx, rIdx, d)
                        * grad_u_srfc<_>(lIdx, rIdx, d, dir);
    }
    return gammaM1() * (grad_rho_E_srfc<_>(lIdx, rIdx, dir)
                        - 0.5 * grad_rho_umag2);
  }

  template<class F>
  auto at_surface(const CellIdx lIdx, const CellIdx rIdx, F&& f) const noexcept {
    return b_()->at_surface(lIdx, rIdx, std::forward<F>(f));
  }


  /// \brief Computes the viscous flux using a three/five-point
  /// central-difference scheme
  ///
  /// $\mathbf{H}_\mathrm{viscous}(\mathbf{Q}) = (0, \mathbf{\tau},
  /// \mathbf{\tau} \mathbf{u} + \mathbf{q})^T$
  template<class _> inline NumA<nvars> compute_num_flux_viscous_
  (const CellIdx lIdx, const CellIdx rIdx, const SInd dir, const Num, const Num,
   flux::viscous_::three_point) const noexcept {
    NumA<nvars> f_v = NumA<nvars>::Zero();

    /// Compute the viscosity and the heat flux:
    const Num mu_srfc = mu(temperature(gamma(), at_surface(lIdx, rIdx, [&](CellIdx c){ return p<_>(c); }),
                                       at_surface(lIdx, rIdx, [&](CellIdx c){ return rho<_>(c); })));
    const Num gradP_srfc = grad_p_srfc<_>(lIdx, rIdx, dir);
    const Num gradRho_srfc = grad_rho_srfc<_>(lIdx, rIdx, dir);

    const Num q = heat_flux(rho_srfc<_>(lIdx, rIdx), gradRho_srfc, mu_srfc,
                            p_srfc<_>(lIdx, rIdx), gradP_srfc);

    /// Compute the shear stress
    auto slope_u = [&, lIdx, rIdx](const SInd c, const SInd slopeDir) {
      return grad_u_srfc<_>(lIdx, rIdx, c, slopeDir);
    };
    NumAM<nd, nd> tau = shear_stress(mu_srfc, slope_u);
    NumA<nd> tau_surface = tau.col(dir);

    /// Viscous flux:
    for (SInd d = 0; d < nd; ++d) {
      f_v(V::rho_u(d)) = 1. / quantities.Re0() * tau_surface(d);
      f_v(V::rho_E()) += u_srfc<_>(lIdx, rIdx, d) * tau_surface(d);
    }
    f_v(V::rho_E()) = 1. / quantities.Re0() * (f_v(V::rho_E()) + q);

    // DBGV((quantities.Re0())(rho_srfc)(p_srfc)(mu_srfc)(gradP_srfc)\
    //     (gradRho_srfc)(q)(shear_stress(mu_srfc, slope_u))(tau_surface)(f_v));
    return f_v;
  }

  ///@}
};

HOM3_FV_PHYSICS_SOLVER_();

}  // namespace cns

////////////////////////////////////////////////////////////////////////////////
}  // namespace fv
}  // namespace solver
}  // namespace hom3
////////////////////////////////////////////////////////////////////////////////
#undef ENABLE_DBG_
////////////////////////////////////////////////////////////////////////////////
#endif
