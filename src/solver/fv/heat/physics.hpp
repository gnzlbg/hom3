#ifndef HOM3_SOLVERS_FV_HEAT_PHYSICS_HPP_
#define HOM3_SOLVERS_FV_HEAT_PHYSICS_HPP_
////////////////////////////////////////////////////////////////////////////////
/// \file \brief Implements the physics class for the Heat-equation.
////////////////////////////////////////////////////////////////////////////////
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

/// \brief Heat-conduction
namespace heat {

/// \brief Implements physics for the Heat-eqt.
template<SInd nd_, class NumFlux, class Solver> struct Physics {
  /// \name Solver Interface
  ///@{
  using physics_type = type_tag;
  static const constexpr SInd nd = nd_;
  static String physics_name() noexcept
  { return "Heat" + std::to_string(nd) + "D"; }
  ///@}

  /// \brief Constructor
  ///
  /// Requires property:
  /// - diffusivity
  /// - T_ref
  /// - CFL
  explicit Physics(io::Properties properties) noexcept
    : quantities(properties)
    , cfl_(io::read<Num>(properties, "CFL"))
  {}

  /// \name Variable Access
  ///@{

  /// Alias for data Indices
  using V = Indices<nd>;

  static constexpr SInd nvars = V::nvars;

  /// \brief Dimensionless temperature
  /// ($T = \overline{T}/\overline{T}_\mathrm{ref} \; [-] $) at cell \p cIdx
  template<class _> inline Num& T(const CellIdx cIdx)       noexcept
  { return b_()->template Q(_(), cIdx, V::T()); }
  template<class _> inline Num  T(const CellIdx cIdx) const noexcept
  { return b_()->template Q(_(), cIdx, V::T()); }
  template<class _> inline Num& T(_, const CellIdx cIdx)       noexcept
  { return b_()->template Q(_(), cIdx, V::T()); }
  template<class _> inline Num  T(_, const CellIdx cIdx) const noexcept
  { return b_()->template Q(_(), cIdx, V::T()); }

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
  ///  dt_i = \frac{ cfl * h_i^2 }{ 2 * mu_i } but mu_i = 1!
  template<class U> inline Num compute_dt(const CellIdx cIdx) const noexcept {
    const Num h = b_()->cells().length(cIdx);
    return cfl_ * 0.5 * std::pow(h, 2.);
  }
  ///@}

  template<class Output> void physics_output(Output&&) const noexcept {}

  template<class _>
  inline void check_variables(const CellIdx cIdx) const noexcept {
    if (b_()->is_ghost_cell(cIdx)) { return; }
    ASSERT(T<_>(cIdx) > 0, "negative temperature (=" << T<_>(cIdx)
                           << ") in cell: " << cIdx << "!\n");
  }

 private:
  /// CRTP:
        Solver* b_()       noexcept { return static_cast<Solver*>(this); }
  const Solver* b_() const noexcept { return static_cast<const Solver*>(this); }

  /// Physical quantities
  const Quantities quantities;

  /// CFL number
  const Num cfl_;

  /// Numerical fluxes implementation
  ///@{

  /// \name Three-point stencil
  ///@{
  template<class _> inline NumA<nvars> compute_num_flux_
  (const CellIdx lIdx, const CellIdx rIdx, const SInd, const Num dx,
    const Num, flux::three_point) const noexcept {
    NumA<nvars> tmp;
    tmp(0) = 1 / dx * (T<_>(lIdx) - T<_>(rIdx));
    return tmp;
  }

  ///@}
  ///@}
};

HOM3_FV_PHYSICS_SOLVER_();

}  // namespace heat

////////////////////////////////////////////////////////////////////////////////
}  // namespace fv
}  // namespace solver
}  // namespace hom3
////////////////////////////////////////////////////////////////////////////////
#undef ENABLE_DBG_
////////////////////////////////////////////////////////////////////////////////
#endif
