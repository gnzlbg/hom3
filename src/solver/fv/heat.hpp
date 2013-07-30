#ifndef HOM3_SOLVERS_FV_HEAT_HPP_
#define HOM3_SOLVERS_FV_HEAT_HPP_
////////////////////////////////////////////////////////////////////////////////
#include "solver.hpp"
/// Options:
#define ENABLE_DBG_ 0
#include "../../misc/dbg.hpp"
////////////////////////////////////////////////////////////////////////////////

namespace solver { namespace fv { namespace heat {

/// \brief Indices to access conservative/primitive variables
template<SInd nd_> struct Indices {
  static constexpr SInd nd = nd_;
  static constexpr SInd nvars = 1;
  static inline constexpr SInd T() { return 0; }
  static inline std::string cv_names(const SInd i) {
    if(i == 0) {
      return std::string("T");
    } else {
      TERMINATE("unknown variable index: " + std::to_string(i));
    }
  }
  static inline std::string pv_names(const SInd i) { return cv_names(i); }
};

struct type_tag {};

namespace flux {
struct three_point {};
}

template<SInd nd, class Solver> struct Physics {

  /// \name Solver Interface
  ///@{

  using type = type_tag;
  static std::string physics_name() { return std::string("Heat") + std::to_string(nd) + "D"; }

  /// \brief Solver type enum (required by solver::Interface)
  static constexpr solver::type type_id() { return solver::type::fv_heat; }

  ///@}

  /// \brief Constructor: requires property conductivity
  Physics(io::Properties properties) : mu_(io::read<Num>(properties,"conductivity")) {}

  /// \name Variable Access
  ///@{

  /// Alias for data Indices
  using V = Indices<nd>;

  static constexpr SInd nvars = V::nvars;

  /// \brief Temperature at cell \p lId
  template<class U> inline Num& T(const Ind lId)       { return Q_<U>(lId,V::T()); }
  template<class U> inline Num  T(const Ind lId) const { return Q_<U>(lId,V::T()); }

  inline Num conductivity() const { return mu_; }

  ///@}

  /// \name Numerical functions
  ///@{

  /// \brief Computes the numerical flux at the interface between \p leftId and
  /// \p rightId
  template<class U, class Flux>
  inline NumA<nvars> compute_num_flux(const Ind leftId, const Ind rightId,
                                      const SInd d, const Num dx, const Num dt) const {
    return compute_num_flux_<U>(leftId,rightId,d,dx,dt,Flux());
  }

  /// \brief computes dt at cell \p lId
  ///  dt_i = \frac{ cfl * h_i^2 }{ 2 * mu_i }
  template<class U> inline Num compute_dt(const Ind lId) const {
    const Num h = b_()->cell_length(lId);
    return b_()->cfl_ * std::pow(h,2) / (2 * conductivity());
  }
  ///@}

  template<class Output> void physics_output(Output&&) const {}

  template<class U> inline void check_variables(const Ind lId) const {
    if(b_()->is_ghost_cell(lId)){ return; }
    ASSERT(T<U>(lId) > 0,   "negative temperature (=" << T<U>(lId) << ") in cell: "   << lId << "!\n");
  }

 private:

  /// CRTP:
        Solver* b_()       { return static_cast<Solver*>(this); }
  const Solver* b_() const { return static_cast<const Solver*>(this); }

  /// Remove template ugliness
  // Quiz fact: template keyword is required because C++ is unparseable -.-
  template<class T> inline Num& Q_(const Ind lId, const SInd varId)       {
    return b_()->template Q<T>(lId,varId);
  }
  template<class T> inline Num  Q_(const Ind lId, const SInd varId) const {
    return b_()->template Q<T>(lId,varId);
  }

  /// Heat conductivity
  const Num mu_;

  /// Numerical fluxes implementation
  ///@{

  /// \name Three-point stencil
  ///@{
  template<class U> inline NumA<nvars> compute_num_flux_
  (const Ind lId, const Ind rId, const SInd, const Num dx,
   const Num, flux::three_point) const {
    NumA<nvars> tmp;
    tmp(0) = conductivity() / dx * ( T<U>(lId) - T<U>(rId) );
    return tmp;
  }

  ///@}
  ///@}
};

template<SInd nd, class NumFlux /* = flux::three_point*/>
using Solver = solver::fv::Solver<nd,Physics,NumFlux>;


namespace bc {

template<class Solver> struct Dirichlet : boundary::Condition {
  template<class F> Dirichlet(Solver& s_, F&& val) : value(val), s(s_) {}

  void apply(AnyRange<Ind> ghost_cells) const override {
    for(auto ghostCellId : ghost_cells) {
      Ind bndryCellId;
      std::tie(bndryCellId,std::ignore) = s.boundary_cell_id(ghostCellId);
      s.cells().vars(ghostCellId,0) = 2.0 * value(bndryCellId,0) - s.cells().vars(bndryCellId,0);
      s.cells().rhs(ghostCellId,0)  = 2.0 * value(bndryCellId,0) - s.cells().rhs (bndryCellId,0);
    }
  }

  std::function<Num(Ind,SInd)> value;
  Solver& s;
};

template<class Solver> struct Neumann : boundary::Condition {
  template<class F> Neumann(Solver& s_, F&& val) : value(val), s(s_) {}

  void apply(AnyRange<Ind> ghost_cells) const override {
    for(auto ghostCellId : ghost_cells) {
      Ind bndryCellId;
      Num dx = s.cells().length(ghostCellId);
      std::tie(bndryCellId,std::ignore) = s.boundary_cell_id(ghostCellId);
      s.cells().vars(ghostCellId,0) = s.cells().vars(bndryCellId,0) - dx * value(bndryCellId,0);
      s.cells().rhs(ghostCellId,0)  = s.cells().rhs (bndryCellId,0) - dx * value(bndryCellId,0);
    }
  }

  std::function<Num(Ind,SInd)> value;
  Solver& s;
};

} // bc namespace

}}} // solver::fv::heat namespace

////////////////////////////////////////////////////////////////////////////////
#undef ENABLE_DBG_
////////////////////////////////////////////////////////////////////////////////
#endif
