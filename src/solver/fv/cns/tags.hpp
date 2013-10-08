#ifndef HOM3_SOLVERS_FV_CNS_TAGS_HPP_
#define HOM3_SOLVERS_FV_CNS_TAGS_HPP_
////////////////////////////////////////////////////////////////////////////////
/// \file \brief Defines the Navier-Stokes equations' tags.
////////////////////////////////////////////////////////////////////////////////
#include "solver/fv/euler/tags.hpp"
////////////////////////////////////////////////////////////////////////////////
namespace hom3 { namespace solver { namespace fv { namespace cns {
////////////////////////////////////////////////////////////////////////////////

/// \brief Navier-Stokes physic's tag
struct type_tag {};

/// \brief Navier Stokes numerical-flux tags
namespace flux {
/// \brief Navier Stokes inviscid flux tags
namespace inviscid_ { using namespace euler::flux; }
/// \brief Navier Stokes viscous flux tags
namespace viscous_ { struct three_point {}; }

/// \brief Standard numerical flux for the Navier-Stokes physics
struct standard : inviscid_::ausm, viscous_::three_point {};
}  // namespace flux

////////////////////////////////////////////////////////////////////////////////
}  // namespace cns
}  // namespace fv
}  // namespace solver
}  // namespace hom3
////////////////////////////////////////////////////////////////////////////////
#endif
