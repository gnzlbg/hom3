#ifndef HOM3_SOLVERS_FV_EULER_PHYSICS_HPP_
#define HOM3_SOLVERS_FV_EULER_PHYSICS_HPP_
////////////////////////////////////////////////////////////////////////////////
/// \file \brief Implements the physics class for the Euler-equations.
////////////////////////////////////////////////////////////////////////////////
#include <limits>
#include <string>
#include <algorithm>
#include "solver/fv/solver.hpp"
#include "indices.hpp"
#include "tags.hpp"
#include "quantities.hpp"
/// Options:
#define ENABLE_DBG_ 0
#include "misc/dbg.hpp"
////////////////////////////////////////////////////////////////////////////////
namespace hom3 { namespace solver { namespace fv {
////////////////////////////////////////////////////////////////////////////////

/// \brief Euler-equations
namespace euler {

/// \brief Implements physics for the Euler-eqts. of compressible
/// gas dynamics.
template<SInd nd_, class NumFlux, class Solver> struct Physics {
  /// \name Solver interface
  ///@{
  using physics_type = type_tag;
  static const constexpr SInd nd = nd_;
  static std::string physics_name() noexcept {
    return "Euler" + std::to_string(nd) + "D";
  }
  ///@}

  /// \brief Constructor:
  ///
  /// Requires property:
  /// - CFL
  explicit Physics(io::Properties properties) noexcept
    : quantities(properties)
    , cfl_(io::read<Num>(properties, "CFL"))
  {}

  /// \name Variable Access
  ///@{

  /// Alias for data Indices
  using V = Indices<nd>;

  /// #of variables: nd + 2 = nd (u_vector) + 1 (rho) + 1 (E)
  static constexpr SInd nvars = V::nvars;

  /// \brief Density of cell \p cIdx
  template<class _> inline Num& rho(const CellIdx cIdx)       noexcept
  { return b_()->Q(_(), cIdx, V::rho()); }
  template<class _> inline Num  rho(const CellIdx cIdx) const noexcept
  { return b_()->Q(_(), cIdx, V::rho()); }

  /// \brief \p d-th momentum component of cell \p cIdx
  template<class _>
  inline Num& rho_u(const CellIdx cIdx, const SInd d)       noexcept
  { return b_()->Q(_(), cIdx, V::rho_u(d)); }
  template<class _>
  inline Num  rho_u(const CellIdx cIdx, const SInd d) const noexcept
  { return b_()->Q(_(), cIdx, V::rho_u(d)); }

  template<class _>
  inline NumA<nd> rho_u(const CellIdx cIdx) const noexcept
  { return b_()->Q(_(), cIdx).template head<nd>(); }

  template<class _>
  inline Eigen::Block<NumM<nd>, 1, nd> rho_u(const CellIdx cIdx) noexcept
  { return b_()->Q(_(), cIdx).template head<nd>(); }

  /// \brief Energy of cell \p cIdx
  template<class _> inline Num& rho_E(const CellIdx cIdx)       noexcept
  { return b_()->Q(_(), cIdx, V::rho_E()); }
  template<class _> inline Num  rho_E(const CellIdx cIdx) const noexcept
  { return b_()->Q(_(), cIdx, V::rho_E()); }

  /// \brief \p d-th velocity component of cell \p cIdx
  template<class _>
  inline Num u(const CellIdx cIdx, const SInd d) const noexcept
  { return rho_u<_>(cIdx, d) / rho<_>(cIdx); }

  /// \brief Velocity vector of cell \p cIdx
  template<class _>
  inline NumA<nd> u(const CellIdx cIdx) const noexcept
  { return rho_u<_>(cIdx).array() / rho<_>(cIdx); }

  /// \brief Pressure of cell \p cIdx p = (\gamma - 1) \rho (E - ||u||_2^2 / 2)
  template<class _> inline Num p(const CellIdx cIdx) const noexcept {
    return gammaM1() * (rho_E<_>(cIdx) - 0.5 * rho<_>(cIdx) * u_mag2<_>(cIdx));
  }

  /// \brief \p d-th component of the Euler flux of cell \p cIdx
  template<class _>
  inline NumA<nvars> flux(const CellIdx cIdx, const SInd d) const noexcept {
    NumA<nvars> f = NumA<nvars>::Zero();

    f(V::rho()) = rho_u<_>(cIdx, d);  // f(rho)= rho u_d

    for (SInd i = 0; i < nd; ++i) {  // rho u_i u_d + kd_di * p
      f(V::rho_u(i)) = rho_u<_>(cIdx, i) * u<_>(cIdx, d);
    }
    f(V::rho_u(d)) += p<_>(cIdx);

    // u_d * (rhoE + p):
    f(V::rho_E()) = u<_>(cIdx, d) * (rho_E<_>(cIdx) + p<_>(cIdx));

    DBGV((cIdx)(d)(f(V::rho()))(f(V::rho_u(0)))(f(V::rho_u(1)))(f(V::rho_E())));
    return f;
  }

  /// \brief Velocity magnitude squared: u^2 + v^2 + w^2 (in 3D)
  template<class _> inline Num u_mag2(const CellIdx cIdx) const noexcept
  { return u<_>(cIdx).squaredNorm(); }

  /// \brief Velocity magnitude: ||u|| = \sqrt{u^2 + v^2 + w^2} (in 3D)
  template<class _> inline Num u_mag(const CellIdx cIdx) const noexcept
  { return u<_>(cIdx).norm(); }

  /// \brief Speed of sound a = \sqrt{\gamma p / \rho}
  template<class _> inline Num a(const CellIdx cIdx) const noexcept
  { return std::sqrt(gamma() * p<_>(cIdx) / rho<_>(cIdx)); }

  /// \brief Mach number at cell \p cIdx (always positive!)
  template<class _> inline Num M(const CellIdx cIdx) const noexcept
  { return u_mag<_>(cIdx) / a<_>(cIdx); }

  /// \brief Specific total Energy $E$ at cell \p cIdx
  template<class _> inline Num E(const CellIdx cIdx) const noexcept
  { return rho_E<_>(cIdx) / rho<_>(cIdx); }

  /// \brief Internal energy density $\rho e$ at cell \p cIdx
  template<class _> inline Num rho_e(const CellIdx cIdx) const noexcept
  { return p<_>(cIdx) / gammaM1(); }

  /// \brief Specific internal energy $e$ at cell \p cIdx
  template<class _> inline Num e(const CellIdx cIdx) const noexcept
  { return rho_e<_>(cIdx) / rho<_>(cIdx); }

  /// \brief Temperature $T$ at cell \p cIdx
  template<class _> inline Num T(const CellIdx cIdx) const noexcept
  { return gamma() * p<_>(cIdx) / rho<_>(cIdx); }

  /// \brief Temperature $T$ at cell \p cIdx
  template<class _> inline Num T(_, const CellIdx cIdx) const noexcept
  { return T<_>(cIdx); }

  inline Num gamma()   const noexcept { return quantities.gamma(); }
  inline Num gammaM1() const noexcept { return quantities.gammaM1(); }

  ///@}

  /// \name Numerical functions
  ///@{

  /// \brief Computes the numerical flux at the interface between \p lIdx and
  /// \p rIdx
  template<class _>
  inline NumA<nvars> compute_num_flux
  (const CellIdx lIdx, const CellIdx rIdx, const SInd d, const Num dx,
   const Num dt) const noexcept
  { return compute_num_flux_<_>(lIdx, rIdx, d, dx, dt, NumFlux()); }

  /// \brief computes dt at cell \p cIdx
  /// \min_{u_i \in \mathbf{u}} ( \frac{C * h_{cell}}{u_i + a} )
  template<class _> inline Num compute_dt(const CellIdx cIdx) const noexcept {
    const Num h = b_()->cells().length(cIdx);
    const Num a_ = a<_>(cIdx);

    /// Compute maximum wave speed: S = max(|u_d| + a)
    const auto S_max
        = boost::accumulate(u<_>(cIdx), std::numeric_limits<Num>::lowest(),
         [a_](Num s, Num u) { return std::max(s, std::abs(u) + a_); });

    const Num dt = cfl_ * h / S_max;

    // Debug dt-kernel:
    DBGV((cIdx)(dt)(a_)(cfl_)(S_max)(h));
    return dt;
  }

  /// \brief Computes the source term
  template<class T>
  inline NumA<nvars> compute_source_term(T, const CellIdx) const noexcept {
    NumA<nvars> source = NumA<nvars>::Zero();
    return source;
  }

  ///@}

  template<class Output> void physics_output(Output&& out) const noexcept {
    out << io::stream("a", 1, [&](const Ind cIdx, const SInd) {
      return a<lhs_tag>(CellIdx{cIdx});
    });
    out << io::stream("p", 1, [&](const Ind cIdx, const SInd) {
      return p<lhs_tag>(CellIdx{cIdx});
    });
    out << io::stream("u", nd, [&](const Ind cIdx, const SInd d) {
      return u<lhs_tag>(CellIdx{cIdx}, d);
    });
    out << io::stream("M", 1, [&](const Ind cIdx, const SInd) {
      return M<lhs_tag>(CellIdx{cIdx});
    });
    out << io::stream("e", 1, [&](const Ind cIdx, const SInd) {
      return e<lhs_tag>(CellIdx{cIdx});
    });
    out << io::stream("T", 1, [&](const Ind cIdx, const SInd) {
      return T<lhs_tag>(CellIdx{cIdx});
    });
  }

  /// \brief checks that variables have only valid values (for debugging
  /// pourposes)
  template<class _>
  inline void check_variables(const CellIdx cIdx) const noexcept {
    if (b_()->is_ghost_cell(cIdx)) { return; }
    ASSERT(rho<_>(cIdx) > 0, "negative density (=" << rho<_>(cIdx)
                             << ") in cell: "  << cIdx << "!\n");
    ASSERT(rho_E<_>(cIdx) > 0,   "negative energy (=" << rho_E<_>(cIdx)
                             << ") in cell: "   << cIdx << "!\n");
    ASSERT(p<_>(cIdx) > 0,   "negative pressure (=" << p<_>(cIdx)
                             << ") in cell: " << cIdx << "!\n");
  }

  /// Physical quantities
  const Quantities quantities;

 protected:
  /// CFL number
  const Num cfl_;

  /// CRTP:
        Solver* b_()       noexcept { return static_cast<Solver*>(this); }
  const Solver* b_() const noexcept { return static_cast<const Solver*>(this); }


  /// Numerical fluxes implementation
  ///@{

  /// \name Local-Lax-Friedrichs Flux
  ///@{

  /// \brief Computes the \p d -th component of the Local-Lax_Friedrichs flux at
  /// the surface between the cells \p lIdx and \p rIdx. The distance between
  /// the cell centers is \p dx and the time-step is \p dt
  template<class _> inline NumA<nvars> compute_num_flux_
  (const CellIdx lIdx, const CellIdx rIdx, const SInd d, const Num dx,
   const Num dt, flux::lax_friedrichs) const noexcept {
    return 0.5 * (flux<_>(lIdx, d) + flux<_>(rIdx, d)
                  + dx / dt * (Q(_(), lIdx) - Q(_(), rIdx)));
  }

  ///@}

  /// \name Advection Upstream Splitting Method (Liu-Steffen 1993)
  ///@{

  /// \brief Computes the \p -dth component of the AUSM flux at the surface
  /// between the cells \p lIdx and \p rIdx.
  template<class _> inline NumA<nvars> compute_num_flux_
  (const CellIdx lIdx, const CellIdx rIdx, const SInd d, const Num,
   const Num, flux::ausm) const noexcept {
    const Num ML = u<_>(lIdx, d) / a<_>(lIdx);
    const Num MR = u<_>(rIdx, d) / a<_>(rIdx);

    const Num m_i = m_int_<+1>(ML) + m_int_<-1>(MR);
    const Num p_i = p_int_<+1>(ML) * p<_>(lIdx) + p_int_<-1>(MR) * p<_>(rIdx);

    NumA<nvars> f_i = m_i >= 0 ? theta_<_>(lIdx) : theta_<_>(rIdx);

    f_i *= m_i;
    f_i(V::rho_u(d)) += p_i;

    DBGV((lIdx)(rIdx)(d)(ML)(MR)(u_mag<_>(lIdx))(u_mag<_>(rIdx))
         (a<_>(lIdx))(a<_>(rIdx))(m_i)(p_i)(f_i));
    return f_i;
  }

 private:
  /// \brief Computes the interface Mach number \mathcal{M}^{+-}
  template<int sign> static inline constexpr Num m_int_(const Num M) noexcept {
    static_assert(sign == 1 || sign == -1, "invalid sign!");
    const Num s = static_cast<Num>(sign);
    return std::abs(M) > 1 ? 0.5 * (M + s * std::abs(M))
                           : s * 0.25 * std::pow(M + s, 2.);
  }

  /// \brief Computes the interface pressure \mathcal{P}^{+-}
  template<int sign> static inline constexpr Num p_int_(const Num M) noexcept {
    static_assert(sign == 1 || sign == -1, "invalid sign!");
    const Num s = static_cast<Num>(sign);
    return std::abs(M) > 1 ? 0.5 * (M + s * std::abs(M))/M
                           : 0.25 * std::pow(M + s, 2.) * (2. - s * M);
  }

  /// \brief Computes the d-th flux vector component divided by the d-th
  /// velocity (it is multiplied by the interface speed-of-sound in ausm)
  template<class _>
  inline NumA<nvars> theta_(const CellIdx cIdx) const noexcept {
    NumA<nvars> f = NumA<nvars>::Zero();
    f(V::rho()) = rho<_>(cIdx);
    for (SInd i = 0; i < nd; ++i) {
      f(V::rho_u(i)) = rho_u<_>(cIdx, i);
    }
    f(V::rho_E()) = rho_E<_>(cIdx) + p<_>(cIdx);

    f *= a<_>(cIdx);
    return f;
  }

  ///@}

  ///@}
 public:

  inline NumA<nvars> cv(const NumA<nvars> pvs) const noexcept
  { return quantities.template cv<nd>(pvs); }

  inline NumA<nvars> pv(const NumA<nvars> cvs) const noexcept
  { return quantities.template pv<nd>(cvs); }

  template<class _> inline NumA<nvars> cv(_, const CellIdx cIdx) const noexcept
  { return b_()->Q(_(), cIdx); }

  template<class _> inline NumA<nvars> pv(_, const CellIdx cIdx) const noexcept
  { return pv(cv(_(), cIdx)); }

};

HOM3_FV_PHYSICS_SOLVER_();

} // namespace euler

////////////////////////////////////////////////////////////////////////////////
} // namespace fv
} // namespace solver
} // namespace hom3
////////////////////////////////////////////////////////////////////////////////
#undef ENABLE_DBG_
////////////////////////////////////////////////////////////////////////////////
#endif
