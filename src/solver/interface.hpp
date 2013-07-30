#ifndef HOM3_SOLVER_INTERFACE_HPP_
#define HOM3_SOLVER_INTERFACE_HPP_
////////////////////////////////////////////////////////////////////////////////
/// \file \brief This file implements the solver interface
////////////////////////////////////////////////////////////////////////////////
#include "../globals.hpp"
////////////////////////////////////////////////////////////////////////////////

/// \brief Solvers
namespace solver {

enum class type {
  fv_heat,
  fv_euler
};

struct Interface {

  template<class S> Interface(const S& s)
      : boundary_condition([&](SInd bcId, Ind localId){ s.boundary_condition(bcId,localId); }),
        solver_id([&](){ return s.solver_id(); }),
        global_ids([&]() { return s.global_ids(); }),
        type_id([&](){ return s.type_id(); })
  {}

  std::function<void(SInd,Ind)> boundary_condition;
  std::function<SInd(void)> solver_id;
  std::function<AnyRange<Ind>()> global_ids;
  std::function<type()> type_id;
};


} // solver namespace
////////////////////////////////////////////////////////////////////////////////
#endif
