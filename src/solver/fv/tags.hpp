#ifndef HOM3_SOLVERS_FV_TAGS_HPP_
#define HOM3_SOLVERS_FV_TAGS_HPP_
////////////////////////////////////////////////////////////////////////////////
/// \brief This file collects the finite volume solver tags
////////////////////////////////////////////////////////////////////////////////
namespace hom3 { namespace solver { namespace fv {

/// \name Tags for rhs and lhs variables
///@{
struct rhs_tag {};
struct lhs_tag {};
constexpr auto rhs = rhs_tag{};
constexpr auto lhs = lhs_tag{};
///@}

namespace bc {

}  // namespace bc

namespace time_integration {

/// \name Tags for time integration
///@{
struct euler_forward {};
struct runge_kutta_2 {};
///@}

}  // namespace time_integration

////////////////////////////////////////////////////////////////////////////////
}  // namespace fv
}  // namespace solver
}  // namespace hom3
////////////////////////////////////////////////////////////////////////////////
#endif
