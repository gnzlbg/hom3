#ifndef HOM3_SOLVERS_FV_HEAT_MOVING_BOUNDARY_CONDITIONS_HPP_
#define HOM3_SOLVERS_FV_HEAT_MOVING_BOUNDARY_CONDITIONS_HPP_
////////////////////////////////////////////////////////////////////////////////
/// \file \brief Implements the Heat Equantion's moving boundary conditions.
////////////////////////////////////////////////////////////////////////////////
/// Includes:
#include "solver/fv/boundary_condition.hpp"
#include "interpolation/rbf.hpp"
#include "geometry/algorithms.hpp"
/// Options:
#define ENABLE_DBG_ 0
#include "misc/dbg.hpp"
////////////////////////////////////////////////////////////////////////////////
namespace hom3 { namespace solver { namespace fv { namespace heat {
////////////////////////////////////////////////////////////////////////////////
namespace bc {
////////////////////////////////////////////////////////////////////////////////

/// \brief Moving boundary conditions
namespace mb {

/// \brief Dirichlet boundary condition for the temperature
template<class Solver> struct Dirichlet : fv::bc::Condition<Dirichlet<Solver>> {
  /// \brief Imposes a \p temperature distribution
  template<class F, class G, class V>
  Dirichlet(Solver& solver, G&& geometry, V&& vel, F&& temperature) noexcept
    : T_srfc(temperature), s(solver)
    , g([=](const NumA<Solver::nd> x) { return (*geometry)(x); })
    , velocity(vel) {}

  /// \brief Imposes a constant surface \p temperature
  template<class G, class V>
  Dirichlet(Solver& solver, G&& geometry, V&& vel, const Num temperature) noexcept
    : T_srfc([=](CellIdx) { return temperature; }), s(solver)
    , g([=](const NumA<Solver::nd> x) { return (*geometry)(x); })
    , velocity(vel) {}

  /// \brief Applies the boundary condition to a range of ghost cells
  template<class _>
  void operator()(_, const CellIdx ghostIdx) const noexcept {
    s.T(_(), ghostIdx)
      = this->mb_dirichlet([&](const CellIdx i) { return s.T(_(), i); },
                           s, g, ghostIdx, T_srfc);
  }

  bool is_moving() const noexcept { return true; }

  /// Surface temperature distribution
  const std::function<Num(CellIdx)> T_srfc;
  Solver& s;
  const std::function<Num(NumA<Solver::nd>)> g;
  const std::function<Num(SInd d)> velocity;
};


/// \brief Neumann boundary condition for the heat_flux
template<class Solver> struct Neumann : fv::bc::Condition<Neumann<Solver>> {
  /// \brief Imposes a \p heat_flux distribution
  template<class F, class G, class V>
  Neumann(Solver& solver, G&& geometry, V&& vel, F&& heat_flux) noexcept
    : flux_srfc(heat_flux), s(solver)
    , g([=](const NumA<Solver::nd> x) { return (*geometry)(x); })
    , velocity(vel) {}

  /// \brief Imposes a constant surface \p heat_flux
  template<class G, class V>
  Neumann(Solver& solver, G&& geometry, V&& vel, const Num heat_flux) noexcept
    : flux_srfc([=](CellIdx) { return heat_flux; }), s(solver)
    , g([=](const NumA<Solver::nd> x) { return (*geometry)(x); })
    , velocity(vel) {}

  /// \brief Applies the boundary condition to a range of ghost cells
  template<class _>
  void operator()(_, const CellIdx ghostIdx) const noexcept {
    s.T(_(), ghostIdx)
      = this->mb_neumann([&](const CellIdx i) { return s.T(_(), i); },
                           s, g, ghostIdx, flux_srfc);
  }

  bool is_moving() const noexcept { return true; }

  /// Surface heat_flux distribution
  const std::function<Num(CellIdx)> flux_srfc;
  Solver& s;
  const std::function<Num(NumA<Solver::nd>)> g;
  const std::function<Num(SInd d)> velocity;
};

}  // namespace mb

////////////////////////////////////////////////////////////////////////////////
}  // namespace bc
}  // namespace heat
}  // namespace fv
}  // namespace solver
}  // namespace hom3
////////////////////////////////////////////////////////////////////////////////
#undef ENABLE_DBG_
////////////////////////////////////////////////////////////////////////////////
#endif
