#ifndef HOM3_SOLVERS_FV_CNS_TAGS_HPP_
#define HOM3_SOLVERS_FV_CNS_TAGS_HPP_
////////////////////////////////////////////////////////////////////////////////
/// \file \brief Defines the Navier-Stokes equations' tags.
////////////////////////////////////////////////////////////////////////////////
#include "solver/fv/euler/tags.hpp"
////////////////////////////////////////////////////////////////////////////////
namespace hom3 { namespace solver { namespace fv { namespace cns {
////////////////////////////////////////////////////////////////////////////////

struct type_tag {};

namespace flux {
namespace inviscid_ { using namespace euler::flux; }
namespace viscous_ { struct three_point {}; }
struct standard : inviscid_::ausm, viscous_::three_point {};
}  // namespace flux

////////////////////////////////////////////////////////////////////////////////
}  // namespace cns
}  // namespace fv
}  // namespace solver
}  // namespace hom3
////////////////////////////////////////////////////////////////////////////////
#endif
