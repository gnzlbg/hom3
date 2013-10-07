#ifndef HOM3_SOLVERS_FV_EULER_INITIAL_CONDITIONS_HPP_
#define HOM3_SOLVERS_FV_EULER_INITIAL_CONDITIONS_HPP_
////////////////////////////////////////////////////////////////////////////////
/// \file \brief Implements the Euler-equantion's initial conditions.
////////////////////////////////////////////////////////////////////////////////
#include "quantities.hpp"
////////////////////////////////////////////////////////////////////////////////
namespace hom3 { namespace solver { namespace fv { namespace euler {
////////////////////////////////////////////////////////////////////////////////

/// \brief Euler-equations initial conditions
namespace ic {

template<SInd nd, EnableIf<traits::equal<SInd, nd, 2>> = traits::dummy>
NumA<nd + 2> isentropic_vortex(const NumA<nd> x, const Num t) noexcept {
  using V = Indices<nd>;
  Num beta = 5;
  Num gamma = 1.4;
  NumA<nd> x0 = {-5.0, -5.0};
  Num u = 1;
  Num v = 1;
  Num x_t = x(0) - u*t;
  Num y_t = x(1) - v*t;

  Num r = std::sqrt(std::pow(x_t - x0(0), 2) + std::pow(y_t - x0(1), 2));

  Num eT1mr2 = std::exp(1 - std::pow(r, 2));
  Num eT21mr2 = std::exp(2 * (1 - std::pow(r, 2)));
  NumA<nd + 2> pvars;

  pvars(V::u(0)) = u - beta * eT1mr2 * (y_t - x0(1)) / (2 * math::pi);
  pvars(V::u(1)) = v + beta * eT1mr2 * (x_t - x0(0)) / (2 * math::pi);
  pvars(V::rho()) = std::pow(
      (1 - ((gamma - 1) * std::pow(beta, 2) * eT21mr2 )
       / (16 * gamma * std::pow(math::pi, 2)) ),
      1 / (gamma-1));
  pvars(V::p()) = std::pow(pvars(V::rho()), gamma);
  return cv<nd>(pvars, 1.4);
}

template<SInd nd, EnableIf<traits::equal<SInd, nd, 3>> = traits::dummy>
NumA<nd + 2> isentropic_vortex(const NumA<nd>, const Num) noexcept {
  NumA<nd + 2> tmp = NumA<nd + 2>::Zero();
  TERMINATE("unimplemented!");
  return tmp;
}

/// \brief Creates an initial condition for a shock-tube problem
///
/// Can be customized along 2 directions only (x0, y1).
template<SInd nd>
auto shock_tube(const Num dir, const Num angle, const Num x0,
                const Num rhoL, const Num umagL, const Num pL,
                const Num rhoR, const Num umagR, const Num pR) {
  auto ic = [=](const NumA<nd> x) {
    using V = Indices<nd>;
    const Num alpha = angle * (2 * math::pi) / 360;
    const Num gamma = 1.4;
    NumA<nd> x_rot;
    x_rot(0) = std::cos(alpha) * x(0) - std::sin(alpha) * x(1);
    x_rot(1) = std::sin(alpha) * x(0) + std::cos(alpha) * x(1);
    NumA<V::nvars> pvars;

    NumA<nd> uL = NumA<nd>::Zero();
    uL(0) = umagL * std::cos(alpha);
    uL(1) = umagL * std::sin(alpha);
    NumA<nd> uR = NumA<nd>::Zero();
    uR(0) = umagR * std::cos(alpha);
    uR(1) = umagR * std::sin(alpha);

    const Num pos = math::approx<Num>(alpha, 0) ? x0 : x0 * std::sqrt(2.);

    if (x_rot(dir) < pos) {
      pvars(V::rho()) = rhoL;
      for (SInd d = 0; d < nd; ++d) {
        pvars(V::u(d)) = uL(d);
      }
      pvars(V::p()) = pL;
    } else {
      pvars(V::rho()) = rhoR;
      for (SInd d = 0; d < nd; ++d) {
        pvars(V::u(d)) = uR(d);
      }
      pvars(V::p()) = pR;
    }
    return cv<nd>(pvars, gamma);
  };

  return ic;
}

}  // namespace ic

////////////////////////////////////////////////////////////////////////////////
}  // namespace euler
}  // namespace fv
}  // namespace solver
}  // namespace hom3
////////////////////////////////////////////////////////////////////////////////
#endif
