#ifndef HOM3_SOLVERS_FV_ADVECTION_INDICES_HPP_
#define HOM3_SOLVERS_FV_ADVECTION_INDICES_HPP_
////////////////////////////////////////////////////////////////////////////////
/// \file \brief Implements Advection-Equation indices.
////////////////////////////////////////////////////////////////////////////////
#include "globals.hpp"
////////////////////////////////////////////////////////////////////////////////
namespace hom3 { namespace solver { namespace fv { namespace advection {
////////////////////////////////////////////////////////////////////////////////

/// \brief Indices to access conservative/primitive variables
template<SInd nd_> struct Indices {
  static const constexpr SInd nd = nd_;
  static const constexpr SInd nvars = 1;
  static inline constexpr SInd q() noexcept { return 0; }
  static inline String cv_names(const SInd i) noexcept {
    if (i == 0) {
      return "q";
    } else {
      TERMINATE("unknown variable index: " + std::to_string(i));
    }
  }
  static inline String pv_names(const SInd i) noexcept
  { return cv_names(i); }
};

////////////////////////////////////////////////////////////////////////////////
}  // namespace advection
}  // namespace fv
}  // namespace solver
}  // namespace hom3
////////////////////////////////////////////////////////////////////////////////
#endif
