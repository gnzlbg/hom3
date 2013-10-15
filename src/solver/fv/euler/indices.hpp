#ifndef HOM3_SOLVERS_FV_EULER_INDICES_HPP_
#define HOM3_SOLVERS_FV_EULER_INDICES_HPP_
////////////////////////////////////////////////////////////////////////////////
/// \file \brief Implements Euler-Equations indices.
////////////////////////////////////////////////////////////////////////////////
#include "globals.hpp"
////////////////////////////////////////////////////////////////////////////////
namespace hom3 { namespace solver { namespace fv { namespace euler {
////////////////////////////////////////////////////////////////////////////////

/// \brief Indices to access conservative/primitive variables
template<SInd nd_> struct Indices {
  static constexpr SInd nd = nd_;
  static constexpr SInd nvars = nd + 2;
  static inline constexpr SInd u(const SInd d) noexcept { return d; }
  static inline constexpr SInd rho_u(const SInd d) noexcept { return d; }
  static inline constexpr SInd rho() noexcept { return nd; }
  static inline constexpr SInd p() noexcept { return nd + 1; }
  static inline constexpr SInd rho_E() noexcept { return nd + 1; }
  static inline String cv_names(const SInd i) noexcept {
    if (i < nd) {
      return "rho_u" + std::to_string(i);
    } else if (i == nd) {
      return "rho";
    } else if (i == nd + 1) {
      return "rho_E";
    } else {
      TERMINATE("unknown variable index: " + std::to_string(i));
    }
  }
  static inline String pv_names(const SInd i) noexcept {
    if (i < nd) {
      return "u" + std::to_string(i);
    } else if (i == nd) {
      return "rho";
    } else if (i == nd + 1) {
      return "p";
    } else {
      TERMINATE("unknown variable index"  + std::to_string(i));
    }
  }
};

////////////////////////////////////////////////////////////////////////////////
}  // namespace euler
}  // namespace fv
}  // namespace solver
}  // namespace hom3
////////////////////////////////////////////////////////////////////////////////
#endif
