#ifndef HOM3_SOLVERS_FV_EULER_HPP_
#define HOM3_SOLVERS_FV_EULER_HPP_
////////////////////////////////////////////////////////////////////////////////
#include "solver.hpp"
/// Options:
#define ENABLE_DBG_ 0
#include "../../misc/dbg.hpp"
////////////////////////////////////////////////////////////////////////////////

namespace hom3 { namespace solver { namespace fv {

namespace euler {

/// \brief Indices to access conservative/primitive variables
template<SInd nd_> struct Indices {
  static constexpr SInd nd = nd_;
  static constexpr SInd nvars = nd + 2;
  static inline constexpr SInd u(const SInd d) { return d; }
  static inline constexpr SInd rho_u(const SInd d) { return d; }
  static inline constexpr SInd rho() { return nd; }
  static inline constexpr SInd p() { return nd + 1; }
  static inline constexpr SInd E() { return nd + 1; }
  static inline std::string cv_names(const SInd i) {
    if(i < nd) {
      return std::string("rho_u") + std::to_string(i);
    } else if (i == nd) {
      return std::string("rho");
    } else if (i == nd + 1) {
      return std::string("E");
    } else {
      TERMINATE("unknown variable index: " + std::to_string(i));
    }
  }
  static inline std::string pv_names(const SInd i) {
    if(i < nd) {
      return std::string("u") + std::to_string(i);
    } else if (i == nd) {
      return std::string("rho");
    } else if (i == nd + 1) {
      return std::string("p");
    } else {
      TERMINATE("unknown variable index"  + std::to_string(i));
    }
  }
};

struct type_tag {};


namespace flux {
struct ausm {};
struct lax_friedrichs {};
}

/// \brief Implements physics class for the Euler-eqts. of compressible
/// gas dynamics.
template<SInd nd, class Solver> struct Physics {

  /// \name Solver interface
  ///@{
  using type = type_tag;
  static std::string physics_name() { return std::string("Euler") + std::to_string(nd) + "D"; }

  /// \brief Solver type enum (required by solver::Interface)
  static constexpr solver::type type_id() { return solver::type::fv_euler; }

  ///@}

  /// Constructor: requires property gamma
  Physics(io::Properties properties)
      : gamma_(io::read<Num>(properties,"gamma")),
        gammaM1_(gamma_ - 1.0)
  {}

  /// \name Variable Access
  ///@{

  /// Alias for data Indices
  using V = Indices<nd>;

  /// #of variables: nd + 2 = nd (u_vector) + 1 (rho) + 1 (E)
  static constexpr SInd nvars = V::nvars;

  /// \brief Density of cell \p lId
  template<class T> inline Num& rho(const CellIdx lId)       {
    return Q_<T>(lId,V::rho());
  }
  template<class T> inline Num  rho(const CellIdx lId) const {
    return Q_<T>(lId,V::rho());
  }

  /// \brief \p d-th momentum component of cell \p lId
  template<class T> inline Num& rho_u(const CellIdx lId, const SInd d)       {
    return Q_<T>(lId,V::rho_u(d));
  }
  template<class T> inline Num  rho_u(const CellIdx lId, const SInd d) const {
    return Q_<T>(lId,V::rho_u(d));
  }

  /// \brief Energy of cell \p lId
  template<class T> inline Num& E(const CellIdx lId)       {
    return Q_<T>(lId,V::E());
  }
  template<class T> inline Num  E(const CellIdx lId) const {
    return Q_<T>(lId,V::E());
  }

  /// \brief \p d-th velocity component of cell \p lId
  template<class T> inline Num u(const CellIdx lId, const SInd d) const {
    return rho_u<T>(lId,d) / rho<T>(lId);
  }

  /// \brief Pressure of cell \p lId p = (\gamma - 1) * (E - \frac{\rho ||u||^2}{2})
  template<class T> inline Num p(const CellIdx lId) const {
    return gammaM1_ * (E<T>(lId) - 0.5 * rho<T>(lId) * u_mag2<T>(lId));
  }

  /// \brief \p d-th component of the Euler flux of cell \p lId
  template<class T> inline NumA<nvars> flux(const CellIdx lId, const SInd d) const {
    NumA<nvars> f = NumA<nvars>::Zero();

    f(V::rho()) = rho_u<T>(lId, d); // f(rho)= rho u_d

    for(SInd i = 0; i < nd; ++i) {
      f(V::rho_u(i)) = rho_u<T>(lId, i) * u<T>(lId, d); // rho u_i u_d + kd_di * p
    }
    f(V::rho_u(d)) += p<T>(lId);

    f(V::E()) = u<T>(lId,d)*(E<T>(lId) + p<T>(lId)); // u_d *(E+p)

    DBGV((lId)(d)(f(V::rho()))(f(V::rho_u(0)))(f(V::rho_u(1)))(f(V::E())));
    return f;
  }

  /// \brief Velocity magnitude squared: u^2 + v^2 + w^2 (in 3D)
  template<class T> inline Num u_mag2(const CellIdx cId) const {
    Num uMag2 = 0;
    for(SInd d = 0; d < nd; ++d) {
      uMag2 += u<T>(cId,d) * u<T>(cId,d);
    }
    return uMag2;
  }

  /// \brief Velocity magnitude: ||u|| = \sqrt{u^2 + v^2 + w^2} (in 3D)
  template<class T> inline Num u_mag(const CellIdx cId) const {
    return std::sqrt(u_mag2<T>(cId));
  }

  /// \brief Speed of sound a = \sqrt{\gamma p / \rho}
  template<class T> inline Num a(const CellIdx cId) const {
    return std::sqrt(gamma() * p<T>(cId) / rho<T>(cId));
  }

  /// \brief Mach number at cell \p cId (always positive!)
  template<class T> inline Num M(const CellIdx cId) const {
    return u_mag<T>(cId) / a<T>(cId);
  }

  inline Num gamma() const { return gamma_; }

  ///@}

  /// \name Numerical functions
  ///@{

  /// \brief Computes the numerical flux at the interface between \p leftId and
  /// \p rightId
  template<class T, class Flux>
  inline NumA<nvars> compute_num_flux(const CellIdx leftId, const CellIdx rightId,
                                      const SInd d, const Num dx, const Num dt) const {
    return compute_num_flux_<T>(leftId,rightId,d,dx,dt,Flux());
  }

  /// \brief computes dt at cell \p lId
  /// \min_{u_i \in \mathbf{u}} ( \frac{C * h_{cell}}{u_i + a} )
  template<class T> inline Num compute_dt(const CellIdx lId) const {
    const Num h = b_()->cells().length(lId);
    const Num a_ = a<T>(lId);

    /// \warning assumming dx is constant...
    Num v = std::numeric_limits<Num>::max();
    for(SInd d = 0; d < nd; ++d) {
      v = std::min(v,std::abs(u<T>(lId,d)) + a_);
    }
    const Num dt = b_()->cfl_ * h / v;

    // Debug dt-kernel:
    DBGV((lId)(dt)(a_)(b_()->cfl_)(v)(h));
    return dt;
  }

  ///@}

  template<class Output> void physics_output(Output&& out) const{
    out << io::stream("a",1,[&](const Ind cId, const SInd){ return a<lhs_tag>(CellIdx{cId}); });
    out << io::stream("p",1,[&](const Ind cId, const SInd){ return p<lhs_tag>(CellIdx{cId}); });
    out << io::stream("u",nd,[&](const Ind cId, const SInd d){ return u<lhs_tag>(CellIdx{cId},d); });
    out << io::stream("M",1,[&](const Ind cId, const SInd){ return M<lhs_tag>(CellIdx{cId}); });
  }

  /// \brief checks that variables have only valid values (for debugging pourposes)
  template<class T> inline void check_variables(const CellIdx lId) const {
    if(b_()->is_ghost_cell(lId)){ return; }
    ASSERT(rho<T>(lId) > 0, "negative density (=" << rho<T>(lId) << ") in cell: "  << lId << "!\n");
    ASSERT(E<T>(lId) > 0,   "negative energy (=" << E<T>(lId) << ") in cell: "   << lId << "!\n");
    ASSERT(p<T>(lId) > 0,   "negative pressure (=" << p<T>(lId) << ") in cell: " << lId << "!\n");
  }

 private:

  /// CRTP:
        Solver* b_()       { return static_cast<Solver*>(this); }
  const Solver* b_() const { return static_cast<const Solver*>(this); }

  /// Remove template ugliness
  // Quiz fact: template keyword is required because C++ is unparseable -.-
  template<class T> inline Num& Q_(const CellIdx lId, const SInd varId)       {
    return b_()->template Q<T>(lId,varId);
  }
  template<class T> inline Num  Q_(const CellIdx lId, const SInd varId) const {
    return b_()->template Q<T>(lId,varId);
  }
  template<class T> inline Eigen::Block<NumM<nvars>,1,nvars> Q_(const CellIdx lId) {
    return b_()->template Q<T>(lId);
  }
  template<class T> inline NumA<nvars> Q_(const CellIdx lId) const {
    return b_()->template Q<T>(lId);
  }

  /// Variables
  const Num gamma_;
  const Num gammaM1_;

  /// Numerical fluxes implementation
  ///@{


  /// \name Local-Lax-Friedrichs Flux
  ///@{

  /// \brief Computes the \p d -th component of the Local-Lax_Friedrichs flux at
  /// the surface between the cells \p lId and \p rId. The distance between the
  /// cell centers is \p dx and the time-step is \p dt
  template<class T> inline NumA<nvars> compute_num_flux_
  (const CellIdx lId, const CellIdx rId, const SInd d, const Num dx,
   const Num dt, flux::lax_friedrichs) const {
    return 0.5*(flux<T>(lId,d) + flux<T>(rId,d) + dx / dt * (Q_<T>(lId) - Q_<T>(rId)));
  }

  ///@}

  /// \name Advection Upstream Splitting Method (Liu-Steffen 1993)
  ///@{

  ///\brief Computes the \p -dth component of the AUSM flux at the surface
  /// between the cells \p lId and \p rId.
  template<class T> inline NumA<nvars> compute_num_flux_
  (const CellIdx lId, const CellIdx rId, const SInd d, const Num, const Num, flux::ausm) const {
    const Num ML = u<T>(lId,d) / a<T>(lId);
    const Num MR = u<T>(rId,d) / a<T>(rId);

    const Num m_i = m_int_<+1>(ML) + m_int_<-1>(MR);
    const Num p_i = p_int_<+1>(ML) * p<T>(lId) + p_int_<-1>(MR) * p<T>(rId);

    NumA<nvars> f_i = m_i >= 0 ? theta_<T>(lId) : theta_<T>(rId);

    f_i *= m_i;
    f_i(V::rho_u(d)) += p_i;

    DBGV((lId)(rId)(d)(ML)(MR)(u_mag<T>(lId))(u_mag<T>(rId))
         (a<T>(lId))(a<T>(rId))(m_i)(p_i)(f_i));
    return f_i;
  }

  /// \brief Computes the interface Mach number \mathcal{M}^{+-}
  template<int sign> inline Num m_int_(const Num M) const {
    static_assert(sign == 1 || sign == -1, "invalid sign!");
    const Num s = static_cast<Num>(sign);
    return std::abs(M) > 1 ? 0.5 * (M + s * std::abs(M))
                           : s * 0.25 * std::pow(M + s,2);
  }

  /// \brief Computes the interface pressure \mathcal{P}^{+-}
  template<int sign> inline Num p_int_(const Num M) const {
    static_assert(sign == 1 || sign == -1, "invalid sign!");
    const Num s = static_cast<Num>(sign);
    return std::abs(M) > 1 ? 0.5 * (M + s * std::abs(M))/M
                           : 0.25 * std::pow(M + s,2) * (2 - s * M);
  }

  /// \brief Computes the d-th flux vector component divided by the d-th
  /// velocity (it is multiplied by the interface speed-of-sound in ausm)
  template<class T> inline NumA<nvars> theta_(const CellIdx cId) const {
    NumA<nvars> f = NumA<nvars>::Zero();
    f(V::rho()) = rho<T>(cId);
    for(SInd i = 0; i < nd; ++i) {
      f(V::rho_u(i)) = rho_u<T>(cId,i);
    }
    f(V::E()) = E<T>(cId) + p<T>(cId);

    f *= a<T>(cId);
    return f;
  }

  ///@}

  ///@}
};

/// \brief Computes conservative variables from vector of primitive variables
template<SInd nd>
static inline NumA<Indices<nd>::nvars> cv(const NumA<Indices<nd>::nvars> pvars, const Num gamma) {
  using V =  Indices<nd>;
  NumA<V::nvars> cvars = NumA<V::nvars>::Zero();
  cvars(V::rho()) = pvars(V::rho());
  Num velMag2 = 0;
  for(SInd d = 0; d < nd; ++d) {
    cvars(V::rho_u(d)) = pvars(V::rho()) * pvars(V::u(d));
    velMag2 += std::pow(pvars(V::u(d)),2);
  }
  cvars(V::E()) = pvars(V::p()) / (gamma - 1.0)
                  + 0.5 * pvars(V::rho()) * velMag2;
  return cvars;
}



template<SInd nd, class NumFlux /* = flux::ausm*/>
using Solver = solver::fv::Solver<nd,Physics,NumFlux>;

template<SInd nd, traits::EnableIf<traits::equal<SInd,nd,2>> = traits::dummy>
NumA<nd + 2> isentropic_vortex(const NumA<nd> x, const Num t) {
  using V = Indices<nd>;
  Num beta = 5;
  Num gamma = 1.4;
  NumA<nd> x0 = {-5.0, -5.0};
  Num u = 1;
  Num v = 1;
  Num x_t = x(0) - u*t;
  Num y_t = x(1) - v*t;

  Num r = std::sqrt( std::pow(x_t - x0(0),2) + std::pow(y_t - x0(1),2) );

  Num eT1mr2 = std::exp(1 - std::pow(r,2));
  Num eT21mr2 = std::exp(2 * (1 - std::pow(r,2)));
  NumA<nd + 2> pvars;

  pvars(V::u(0)) = u - beta * eT1mr2 * (y_t - x0(1)) / (2 * math::pi);
  pvars(V::u(1)) = v + beta * eT1mr2 * (x_t - x0(0)) / (2 * math::pi);
  pvars(V::rho()) = std::pow(
      (1 - ((gamma - 1) * std::pow(beta,2) * eT21mr2 )/ (16 * gamma * std::pow(math::pi,2)) ),
      1 / (gamma-1));
  pvars(V::p()) = std::pow(pvars(V::rho()),gamma);
  return cv<nd>(pvars,1.4);
}

template<SInd nd, traits::EnableIf<traits::equal<SInd,nd,3>> = traits::dummy>
NumA<nd + 2> isentropic_vortex(const NumA<nd>, const Num) {
  NumA<nd+2> tmp = NumA<nd+2>::Zero();
  TERMINATE("unimplemented!");
  return tmp;
}


namespace bc {

template<class Solver> struct Neumann : boundary::Condition {
  using F = std::function<Num(CellIdx,SInd)>;
  Neumann(Solver& s_, F&& val = [](CellIdx,SInd){return 0.0; }) : value(val), s(s_) {}
  void apply(AnyRange<CellIdx> ghost_cells) const override {
    for(auto ghostCellId : ghost_cells) {
      CellIdx bndryCellId;
      std::tie(bndryCellId,std::ignore) = s.boundary_cell_id(ghostCellId);
      for(SInd v = 0; v < Solver::nvars; ++v) {
        s.cells().lhs(ghostCellId,v) = s.cells().lhs(bndryCellId,v) /* + ...*/;
        s.cells().rhs(ghostCellId,v) = s.cells().rhs(bndryCellId,v) /* + ...*/;
      }
    }
  }

  F value;
  Solver& s;
};

template<class Solver> struct Test : boundary::Condition {
  using V = Indices<Solver::nd>;
  using F = std::function<Num(CellIdx,SInd)>;
  Test(Solver& s_, F&& val) : value(val), s(s_) {}

  void apply(AnyRange<CellIdx> ghost_cells) const override {
    for(auto ghostCellId : ghost_cells) {
       for(SInd v = 0; v < Solver::nvars; ++v) {
         s.cells().lhs(ghostCellId,v) = value(0,v);
       }
    }
  }

  F value;
  Solver& s;
};

template<class Solver> struct IsentropicVortex : boundary::Condition {
  using V = Indices<Solver::nd>;
  IsentropicVortex(Solver& s_) : s(s_) {}

  void apply(AnyRange<CellIdx> ghost_cells) const override {
    for(auto ghostCellId : ghost_cells) {
      CellIdx bndryCellId;
      std::tie(bndryCellId,std::ignore) = s.boundary_cell_id(ghostCellId);

      NumA<V::nd>    x_bc = s.x_cell(bndryCellId);
      NumA<V::nd>    x_gc = s.x_cell(ghostCellId);
      NumA<V::nvars> gcVars = isentropic_vortex<V::nd>(x_gc,s.time());
      NumA<V::nvars> bcVars = isentropic_vortex<V::nd>(x_bc,s.time());
      /// \todo vars = vars
      for(SInd v = 0; v < V::nvars; ++v) {
        s.cells().lhs(ghostCellId,v) = gcVars(v);
        s.cells().rhs(ghostCellId,v) = gcVars(v);
        s.cells().lhs(bndryCellId,v) = bcVars(v);
        s.cells().rhs(bndryCellId,v) = bcVars(v);
      }
    }
  }

  Solver& s;
};


} // bc namespace

}}}} // hom3::solver::fv::euler namespace

////////////////////////////////////////////////////////////////////////////////
#undef ENABLE_DBG_
////////////////////////////////////////////////////////////////////////////////
#endif
