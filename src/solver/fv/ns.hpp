#ifndef HOM3_SOLVERS_FV_CNS_HPP_
#define HOM3_SOLVERS_FV_CNS_HPP_
////////////////////////////////////////////////////////////////////////////////
#include <limits>
#include <string>
#include <algorithm>
#include "solver.hpp"
/// Options:
#define ENABLE_DBG_ 0
#include "../../misc/dbg.hpp"
////////////////////////////////////////////////////////////////////////////////

namespace solver { namespace fv { namespace cns {

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

namespace detail_ {

namespace inviscid {
struct ausm {};
} // inviscid namespace

namespace viscous {
struct three_point {};
} // viscous namespace

} // detail_ namespace

struct standard : detail_::inviscid::ausm, detail_::viscous::three_point {};
}

/// \brief Implements physics class for the Navier-Stokes eqts for compressible
/// flow.
template<SInd nd, class Solver> struct Physics {

  /// \name Solver interface
  ///@{
  using type = type_tag;
  static std::string physics_name() { return std::string("CNS") + std::to_string(nd) + "D"; }

  /// \brief Solver type enum (required by solver::Interface)
  static constexpr solver::type type_id() { return solver::type::fv_cns; }

  ///@}

  /// Constructor: requires property gamma
  Physics(io::Properties properties) :
      /// Read input variables:
      gamma_   (io::read<Num>(properties,"gamma")  ),
      ReInfty_ (io::read<Num>(properties,"ReInfty")),
      MInfty_  (io::read<Num>(properties,"MaInfty")),
      Pr0_     (io::read<Num>(properties,"Pr0")    ),
      Lref_    (io::read<Num>(properties,"Lref")   ),
      Tref_    (io::read<Num>(properties,"Tref")   ),
      S_       (io::read<Num>(properties,"S")      ),
      R_       (io::read<Num>(properties,"R")      ),
      /// Set internal variables
      gammaM1_             (gamma_ - 1.0),
      sOverTref_           (S_/Tref_),
      sOverTrefP1_         (sOverTref_ + 1.0),
      TInfty_              (1 / (1 + 0.5 * gammaM1_ * std::pow(MInfty_,2))),
      rhoInfty_            (std::pow(TInfty, 1 / gammaM1_)),
      pInfty_              (1 / gamma_ * std::pow(TInfty_, gamma_ / gammaM1_)),
      uInfty_              (MInfty_ * std::sqrt(TInfty_)),
      muInfty_             (sutherlands_law(TInfty_)),
      Re0_                 (ReInfty_ * muInfty_ / (rhoInfty_ * MInfty_ * std::sqrt(TInfty_))),
      dt_stagnation2infty_ (uInfty_)
  {}

  /// \name Variable Access
  ///@{

  /// Alias for data Indices
  using V = Indices<nd>;

  /// #of variables: nd + 2 = nd (u_vector) + 1 (rho) + 1 (E)
  static constexpr SInd nvars = V::nvars;

  /// \brief Density of cell \p lId
  template<class T> inline Num& rho(const Ind lId)       {
    return Q_<T>(lId,V::rho());
  }
  template<class T> inline Num  rho(const Ind lId) const {
    return Q_<T>(lId,V::rho());
  }

  /// \brief \p d-th momentum component of cell \p lId
  template<class T> inline Num& rho_u(const Ind lId, const SInd d)       {
    return Q_<T>(lId,V::rho_u(d));
  }
  template<class T> inline Num  rho_u(const Ind lId, const SInd d) const {
    return Q_<T>(lId,V::rho_u(d));
  }

  /// \brief Energy of cell \p lId
  template<class T> inline Num& E(const Ind lId)       {
    return Q_<T>(lId,V::E());
  }
  template<class T> inline Num  E(const Ind lId) const {
    return Q_<T>(lId,V::E());
  }

  /// \brief \p d-th velocity component of cell \p lId
  template<class T> inline Num u(const Ind lId, const SInd d) const {
    return rho_u<T>(lId,d) / rho<T>(lId);
  }

  /// \brief Pressure of cell \p lId p = (\gamma - 1) * (E - \frac{\rho ||u||^2}{2})
  template<class T> inline Num p(const Ind lId) const {
    return gammaM1_ * (E<T>(lId) - 0.5 * rho<T>(lId) * u_mag2<T>(lId));
  }

  /// \brief Velocity magnitude squared: u^2 + v^2 + w^2 (in 3D)
  template<class T> inline Num u_mag2(const Ind cId) const {
    Num uMag2 = 0;
    for(SInd d = 0; d < nd; ++d) {
      uMag2 += u<T>(cId,d) * u<T>(cId,d);
    }
    return uMag2;
  }

  /// \brief Velocity magnitude: ||u|| = \sqrt{u^2 + v^2 + w^2} (in 3D)
  template<class T> inline Num u_mag(const Ind cId) const {
    return std::sqrt(u_mag2<T>(cId));
  }

  /// \brief Speed of sound a = \sqrt{\gamma p / \rho}
  template<class T> inline Num a(const Ind cId) const {
    return std::sqrt(gamma() * p<T>(cId) / rho<T>(cId));
  }

  /// \brief Mach number at cell \p cId (always positive!)
  template<class T> inline Num M(const Ind cId) const {
    return u_mag<T>(cId) / a<T>(cId);
  }

  inline Num gamma() const { return gamma_; }
  inline Num gammaM1() const { return gammaM1_; }


  /// \brief Computes viscosity according to Sutherland's law for the
  /// temperature \param T
  ///
  /// Sutherland's law:
  /// \mu (T) = (T)^{3/2} * \frac{1 + S / T_{\mathrm{ref}}}{T + S / T_{\mathrm{ref}}}
  ///
  /// Here sOverTref_ = S / T_ref , and
  ///      sOverTrefP1_ = 1 + S / T_ref
  ///
  inline Num sutherlands_law(const Num T) const {
    return std::pow(T,1.5) * sOverTrefP1_ / (T + sOverTref_);
  }

  /// \brief Computes the temperature gradient
  ///
  /// \nabla T = \frac{\gamma}{\rho} \left ( \nabla p - \frac{p}{\rho} \nabla \rho \right )
  inline Num grad_T(const Num rho, const Num gradRho,
                    const Num p, const Num gradP) const {
    return gamma() / rho * (gradP - p / rho * gradRho);
  }

  /// \brief Computes the heat flux
  ///
  /// Heat flux:
  /// q = \frac{\mu \gamma}{(\gamma - 1) \mathrm{Pr}_0 \rho}
  /// \left ( \nabla p - \frac{p}{\rho} \nabla \rho \right )
  inline Num heat_flux(const Num rho, const Num gradRho, const Num mu,
                       const Num p, const Num gradP) const {
    return - mu * gamma() / (gammaM1() * Pr0_ * rho) * gradT();
  }

  /// \brief Computes the stress tensor
  ///
  /// Stress-tensor
  /// \tau = \mu (- (\nabla \mathbf{u} + (\nabla \mathbf{u})^{T})
  ///             + \frac{2}{3} (\nabla \cdot \mathbf{u}) \mathbf{I})
  inline NumAM<nd,nd> stress_tensor
  (const Num mu, const NumAM<nd,nd> covD_u, const Num div_u) const {
    return mu * (-(covD_u + covD_u.transpose())
                 + 2/3 * div_u * NumAM<nd,nd>::Identity());
  }

  ///@}

  /// \name Numerical functions
  ///@{

  /// \brief Computes the numerical flux at the interface between \p leftId and
  /// \p rightId
  template<class T, class Flux>
  inline NumA<nvars> compute_num_flux(const Ind leftId, const Ind rightId,
                                      const SInd d, const Num dx, const Num dt) const {
    return compute_num_flux_inviscid_<T>(leftId,rightId,d,dx,dt,Flux())
           + compute_num_flux_viscous_<T>(leftId,rightId,d,dx,dt,Flux());
  }

  /// \brief computes dt at cell \p lId
  /// \min_{u_i \in \mathbf{u}} ( \frac{C * h_{cell}}{u_i + a} )
  template<class T> inline Num compute_dt(const Ind lId) const {
    const Num h = b_()->cell_length(lId);
    const Num a_ = a<T>(lId);

    /// \warning assumming dx is constant...
    Num v = std::numeric_limits<Num>::max();
    for(SInd d = 0; d < nd; ++d) {
      v = std::min(v,std::abs(u<T>(lId,d)) + a_);
    }
    const Num dt = b_()->cfl_ * h / v;

    TERMINATE("viscous time step not implemented");

    // Debug dt-kernel:
    DBGV((lId)(dt)(a_)(b_()->cfl_)(v)(h));
    return dt;
  }

  ///@}

  template<class Output> void physics_output(Output&& out) const{
    out << io::stream("a",1,[&](const Ind cId, const SInd){ return a<lhs_tag>(cId); });
    out << io::stream("p",1,[&](const Ind cId, const SInd){ return p<lhs_tag>(cId); });
    out << io::stream("u",nd,[&](const Ind cId, const SInd d){ return u<lhs_tag>(cId,d); });
    out << io::stream("M",1,[&](const Ind cId, const SInd){ return M<lhs_tag>(cId); });
  }

  /// \brief checks that variables have only valid values (for debugging pourposes)
  template<class T> inline void check_variables(const Ind lId) const {
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
  template<class T> inline Num& Q_(const Ind lId, const SInd varId)       {
    return b_()->template Q<T>(lId,varId);
  }
  template<class T> inline Num  Q_(const Ind lId, const SInd varId) const {
    return b_()->template Q<T>(lId,varId);
  }
  template<class T> inline Eigen::Block<NumM<nvars>,1,nvars> Q_(const Ind lId) {
    return b_()->template Q<T>(lId);
  }
  template<class T> inline NumA<nvars> Q_(const Ind lId) const {
    return b_()->template Q<T>(lId);
  }

  /// \name Variables
  ///@{

  /// \name Input variables
  ///@{

  /// Ratio of specific heats: \gamma = \overline{c}_p / \overline{c}_V
  const Num gamma_;

  /// Free stream Reynolds' number: \mathrm{Re}_\infty = \frac{ \overline{\rho}_\infty \lVert \overline{\mathbf{u}}_\infty \rVert_{L2} \overline{L}_\mathrm{ref} }{ \overline{\mu}_\infty }
  const Num ReInfty_;

  /// Free stream Mach's number: M_\infty = \frac{ \lVert \overline{\mathbf{u}} \rVert_{L2} }{ \overline{a}_\infty }
  const Num MInfty_;

  /// Reference temperature
  const Num Tref_;

  /// Reference length
  const Num Lref_;

  /// Prandtl number: \mathrm{Pr}_{0} = \frac{ \mu_0 c_{p,0} }{ \lambda_{0} }
  const Num Pr0_;

  /// Universal gas constant
  const Num R_;

  ///@}

  /// \name Internal variables
  ///@{

  /// Gamma minus one: \gamma - 1
  const Num gammaM1_;

  /// Free-stream temperature T_\infty =  \left ( 1 + \frac{\gamma - 1}{2} M_\infty^2 \right )^{-1}
  const Num TInfty_;

  /// Free-stream density \rho_\infty = (T_\infty)^{\frac{1}{\gamma - 1}}
  const Num rhoInfty_;

  /// Free-stream pressure p_\infty = \frac{1}{\gamma}(T_\infty)^{\frac{\gamma}{\gamma - 1}}
  const Num pInfty_;

  /// Free-stream velocity u_\infty = M_\infty \sqrt{T_\infty}
  const Num uInfty_;

  /// Free-stream kinematic viscosity \mu_\infty = (T_\infty)^{3}{2} \frac{1 + S/T_{\mathrm{ref}}}{T_\infty + S/T_{\mathrm{ref}}}
  const Num muInfty_;

  /// Stagnation Reynold's number: \mathrm{Re}_0 = \frac{ \mathrm{Re}_\infty \mu_\infty }{ \rho_\infty M_\infty \sqrt{T_\infty} }
  const Num Re0_;

  /// Free-stream time-step conversion factor: factor = u_\infty
  const Num stagnation2infty_dt_;

  ///@}

  ///@}

  /// Numerical fluxes implementation
  ///@{

  /// \name Advection Upstream Splitting Method (Liu-Steffen 1993)
  ///@{

  ///\brief Computes the \p -dth component of the AUSM flux at the surface
  /// between the cells \p lId and \p rId.
  template<class T> inline NumA<nvars> compute_num_flux_inviscid_
  (const Ind lId, const Ind rId, const SInd d, const Num, const Num, flux::ausm) const {
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
  template<class T> inline NumA<nvars> theta_(const Ind cId) const {
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

  /// \brief Computes the viscous flux using a central-difference scheme (three-point stencil)
  template<class T> inline NumA<nvars> compute_num_flux_viscous_
  (const Ind lId, const Ind rId, const SInd d, const Num dx, const Num, flux::three_point) const {
    NumA<nvars> f_v = NumA<nvars>::Zero();

    const NumA<nd> u_srfc = NumA<nd>::Zero();
    for(const SInd d = 0; d < nd; ++d) {
      u_srfc(d) = 0.5 * (u<T>(lId,d) + u<T>(rId,d));
    }
    const Num rho_srfc = 0.5 * (rho<T>(lId) + rho<T>(rId));
    const Num p_srfc = 0.5 * (p<T>(lId) + p<T>(rId));
    const Num T_srfc = gamma() * p_srfc / rho_srfc;
    const Num mu_srfc = sutherlands_law(T_srfc);

    // heat flux:
    const Num gradP_srfc = central_diff(p<T>(lId), p<T>(rId), dx);
    const Num gradRho_srfc = central_diff(rho<T>(lId), rho<T>(rId), dx);
    const Num q = mu_srfc * gamma() / (gammaM1() * Pr0_* rho_srfc)
                  * (gradP_srfc - p_srfc / rho_srfc * gradRho_srfc);

    // velocity gradient
    const NumA<nd> grad_u_srfc = central_diff(u<T>(lId),u<T>(rId),NumA<nd>::Constant(dx));
    
    // stress tensor:
    NumA<nd> tau;
    if(nd == 2) {
      tau(0) = mu_srfc * (4/3 * grad_u_srfc
    } else if (nd == 3) {

    } else {
      static_assert(nd == 2 || nd == 3, "Unsupported dimension!");
    }

    return f_v;
  }

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

template<SInd nd, EnableIf<traits::equal<SInd,nd,2>> = traits::dummy>
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

template<SInd nd, EnableIf<traits::equal<SInd,nd,3>> = traits::dummy>
NumA<nd + 2> isentropic_vortex(const NumA<nd>, const Num) {
  NumA<nd+2> tmp = NumA<nd+2>::Zero();
  TERMINATE("unimplemented!");
  return tmp;
}


namespace bc {

template<class Solver> struct Neumann : boundary::Condition {
  template<class F> Neumann(Solver& s_, F&& val) : value(val), s(s_) {}
  void apply(AnyRange<Ind> ghost_cells) const override {
    for(auto ghostCellId : ghost_cells) {
      Ind bndryCellId;
      std::tie(bndryCellId,std::ignore) = s.boundary_cell_id(ghostCellId);
      for(SInd v = 0; v < Solver::nvars; ++v) {
        s.cells().vars(ghostCellId,v) = s.cells().vars(bndryCellId,v) /* + ...*/;
        s.cells().rhs(ghostCellId,v) = s.cells().rhs(bndryCellId,v) /* + ...*/;
      }
    }
  }

  std::function<Num(Ind,SInd)> value;
  Solver& s;
};

template<class Solver> struct Test : boundary::Condition {
  using V = Indices<Solver::nd>;
  template<class F> Test(Solver& s_, F&& val) : value(val), s(s_) {}

  void apply(AnyRange<Ind> ghost_cells) const override {
    for(auto ghostCellId : ghost_cells) {
       for(SInd v = 0; v < Solver::nvars; ++v) {
         s.cells().vars(ghostCellId,v) = value(0,v);
       }
    }
  }

  std::function<Num(Ind,SInd)> value;
  Solver& s;
};

template<class Solver> struct IsentropicVortex : boundary::Condition {
  using V = Indices<Solver::nd>;
  IsentropicVortex(Solver& s_) : s(s_) {}

  void apply(AnyRange<Ind> ghost_cells) const override {
    for(auto ghostCellId : ghost_cells) {
      Ind bndryCellId;
      std::tie(bndryCellId,std::ignore) = s.boundary_cell_id(ghostCellId);

      NumA<V::nd>    x_bc = s.cells().xc().row(bndryCellId);
      NumA<V::nd>    x_gc = s.cells().xc().row(ghostCellId);
      NumA<V::nvars> gcVars = isentropic_vortex<V::nd>(x_gc,s.time());
      NumA<V::nvars> bcVars = isentropic_vortex<V::nd>(x_bc,s.time());
      for(SInd v = 0; v < V::nvars; ++v) {
        s.cells().vars(ghostCellId,v) = gcVars(v);
        s.cells().rhs(ghostCellId,v) = gcVars(v);
        s.cells().vars(bndryCellId,v) = bcVars(v);
        s.cells().rhs(bndryCellId,v) = bcVars(v);
      }
    }
  }

  Solver& s;
};


} // bc namespace

}}} // solver::fv::euler namespace

////////////////////////////////////////////////////////////////////////////////
#undef ENABLE_DBG_
////////////////////////////////////////////////////////////////////////////////
#endif
