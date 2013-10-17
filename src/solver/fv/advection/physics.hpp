#ifndef HOM3_SOLVERS_FV_ADVECTION_PHYSICS_HPP_
#define HOM3_SOLVERS_FV_ADVECTION_PHYSICS_HPP_
////////////////////////////////////////////////////////////////////////////////
/// \file \brief Implements the physics class for the Advection equation.
////////////////////////////////////////////////////////////////////////////////
#include <algorithm>
#include <limits>
#include "solver/fv/solver.hpp"
#include "indices.hpp"
#include "tags.hpp"
/// Options:
#define ENABLE_DBG_ 0
#include "misc/dbg.hpp"
////////////////////////////////////////////////////////////////////////////////
namespace hom3 { namespace solver { namespace fv {
////////////////////////////////////////////////////////////////////////////////

/// \brief Advection equation
namespace advection {

/// \brief Implements physics for the Advection eqt.
template<SInd nd_, class NumFlux, class Solver> struct Physics {
  /// \name Solver Interface
  ///@{
  using physics_type = type_tag;
  static const constexpr SInd nd = nd_;
  static String physics_name() noexcept
  { return "Advection" + std::to_string(nd) + "D"; }
  ///@}

  /// \brief Constructor
  ///
  /// Requires property:
  /// - CFL
  explicit Physics(io::Properties properties) noexcept
    : cfl_(io::read<Num>(properties, "CFL"))
    , u_(io::read<std::function<NumA<nd>(NumA<nd>)>>(properties, "velocity"))
  {}

  /// \name Variable Access
  ///@{

  /// Alias for data Indices
  using V = Indices<nd>;

  static constexpr SInd nvars = V::nvars;

  /// \brief Access quantity
  template<class _> inline Num& q(const CellIdx cIdx)       noexcept
  { return b_()->template Q(_(), cIdx, V::q()); }
  template<class _> inline Num  q(const CellIdx cIdx) const noexcept
  { return b_()->template Q(_(), cIdx, V::q()); }
  template<class _> inline Num& q(_, const CellIdx cIdx)       noexcept
  { return b_()->template Q(_(), cIdx, V::q()); }
  template<class _> inline Num  q(_, const CellIdx cIdx) const noexcept
  { return b_()->template Q(_(), cIdx, V::q()); }

  ///@}

  /// \name Numerical functions
  ///@{

  /// \brief Computes the numerical flux at the interface between \p leftId and
  /// \p rightId
  template<class U>
  inline NumA<nvars> compute_num_flux
  (const CellIdx lIdx, const CellIdx rIdx, const SInd d, const Num dx,
    const Num dt) const noexcept
  { return compute_num_flux_<U>(lIdx, rIdx, d, dx, dt, NumFlux()); }

  /// \brief Computes the source term
  template<class T>
  inline NumA<nvars> compute_source_term(T, const CellIdx) const noexcept {
    NumA<nvars> source = NumA<nvars>::Zero();
    return source;
  }

  /// \brief computes dt at cell \p cIdx
  ///  dt_i = \frac{ cfl * h_i }{ u }
  template<class U> inline Num compute_dt(const CellIdx cIdx) const noexcept {
    const Num h = b_()->cells().length(cIdx);
    const NumA<nd> xc = b_()->cells().x_center.row(cIdx);
    const auto u_max
      = boost::accumulate(u_(xc), std::numeric_limits<Num>::lowest(),
      [](Num o, Num u) { return std::max(o, std::abs(u)); });
    return cfl_ * h / u_max;
  }
  ///@}

  template<class Output> void physics_output(Output&& out) const noexcept {
    out << io::stream("u", nd, [&](const Ind cIdx, const SInd d) -> Num {
        const NumA<nd> x_c = b_()->cells().x_center.row(CellIdx{cIdx});
        return u_(x_c)(d);
    });
  }

  template<class _>
  inline void check_variables(const CellIdx cIdx) const noexcept {
    if (b_()->is_ghost_cell(cIdx)) { return; }
    ASSERT(q<_>(cIdx) > 0, "negative quantity (=" << q<_>(cIdx)
           << ") in cell: " << cIdx << "!\n");
  }

  /// \brief Advection velocity distribution at \p cIdx
  inline NumA<nd> velocity(CellIdx cIdx) const noexcept
  { return u_(b_()->cells().x_center.row(cIdx)); }

 private:
  /// CRTP:
        Solver* b_()       noexcept { return static_cast<Solver*>(this); }
  const Solver* b_() const noexcept { return static_cast<const Solver*>(this); }

  /// CFL number
  const Num cfl_;
  const std::function<NumA<nd>(NumA<nd>)> u_;

  /// Numerical fluxes implementation
  ///@{

  /// \brief Local Lax-Friedrichs flux
  template<class _> inline NumA<nvars> compute_num_flux_
  (const CellIdx lIdx, const CellIdx rIdx, const SInd d, const Num,
  const Num, flux::local_lax_friedrichs) const noexcept {
    NumA<nvars> tmp;
    const Num qL = q<_>(lIdx), qR = q<_>(rIdx);
    const NumA<nd> uL = velocity(lIdx), uR = velocity(rIdx);
    const Num lambda_max = std::max(std::abs(uL(d)), std::abs(uR(d)));
    tmp(0) = 0.5 * (qL * uL(d) + qR * uR(d) + lambda_max * (qL - qR));
    return tmp;
  }

  /// \brief Upwind flux
  template<class _> inline NumA<nvars> compute_num_flux_
  (const CellIdx lIdx, const CellIdx rIdx, const SInd d, const Num,
  const Num, flux::upwind) const noexcept {
    NumA<nvars> tmp;
    const Num uL = velocity(lIdx)(d),  uR = velocity(rIdx)(d);
    const Num u_srfc = 0.5 * (uL + uR);

    tmp(0) = uL * q<_>(lIdx);
    if (u_srfc < 0) {
      tmp(0) *= -1.;
    }
    return tmp;
  }

  ///@}
  ///@}
};

HOM3_FV_PHYSICS_SOLVER_();

}  // namespace advection

////////////////////////////////////////////////////////////////////////////////
}  // namespace fv
}  // namespace solver
}  // namespace hom3
////////////////////////////////////////////////////////////////////////////////
#undef ENABLE_DBG_
////////////////////////////////////////////////////////////////////////////////
#endif
