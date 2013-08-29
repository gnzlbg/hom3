#ifndef HOM3_SOLVERS_FV_UTILITIES_HPP_
#define HOM3_SOLVERS_FV_UTILITIES_HPP_
////////////////////////////////////////////////////////////////////////////////
/// \brief This file collects some finite volume utilitie functions that are
/// helpful for creating, runing, and managing FV solvers.
////////////////////////////////////////////////////////////////////////////////
#include "globals.hpp"
#include "grid/grid.hpp"
#include "grid/helpers.hpp"
#include "solver.hpp"
////////////////////////////////////////////////////////////////////////////////
namespace hom3 { namespace solver { namespace fv {

/// \brief Add unit cube bcs to solver
template<class Solver, class RootCell, class BCs>
void append_bcs(Solver& solver, RootCell rCell, BCs&& cs) {
  using namespace grid::helpers; using namespace cube; using namespace edge;
  static const SInd nd = Solver::nd;
  auto bcs = make_boundaries<nd>()(solver, rCell, std::forward<BCs>(cs));
  solver.append_bcs(bcs);
}

/// \brief Initializes an object that models the Initializable concept
/// (i.e. provides the initialize() member function)
template<class Initializable>
inline void initialize(Initializable&& initializable) noexcept
{ initializable.initialize(); }
template<class Initializable, class... Initializables>
inline void initialize
(Initializable&& initializable, Initializables&&... initializables) noexcept {
  initialize(std::forward<Initializable>(initializable));
  initialize(std::forward<Initializables>(initializables)...);
}

/// \brief Writes an object that models the Domain concept
/// (i.e. can be written using write_domain(domain)).
template<class Domain> inline void write_domains(Domain& domain) noexcept
{ write_domain(domain); }
template<class Domain, class... Domains> inline void write_domains
(Domain& domain, Domains&... domains) noexcept {
  write_domains(domain);
  write_domains(domains...);
}

template<class Solver>
inline void write_output(const Ind outputInterval, Solver& solver) {
  if(solver.step() % outputInterval == 0) {
    write_domain(solver);
  }
}
template<class Solver, class... Solvers> inline void write_output
(const Ind outputInterval, Solver& solver, Solvers&... solvers) {
  write_output(outputInterval, solver);
  write_output(outputInterval, solvers...);
}

template<class Solver>
bool solver_finished(const Ind maxNoSteps, const Solver& solver) noexcept {
  return solver.time() >= solver.final_time() || solver.step() >= maxNoSteps;
}
template<class Solver, class... Solvers> bool solver_finished
(const Ind maxNoSteps, const Solver& solver, const Solvers&... solvers) noexcept {
  return solver_finished(maxNoSteps, solver)
      || solver_finished(maxNoSteps, solvers...);
}

template<class Solver>
bool solution_diverged(const Solver& solver) noexcept {
  return solver.dt() < 1e-300 || solver.dt() > 1e300;
}

template<class Solver, class... Solvers>
bool solution_diverged(const Solver& solver, const Solvers&... solvers) {
  return solution_diverged(solver) || solution_diverged(solvers...);
}

template<class Solver>
void write_timestep(const Solver& solver) noexcept {
  std::cerr << "Solver: " << solver.domain_name() << " | "
            << "step: " << solver.step() << " | "
            << "time: " << solver.time() << " | "
            << "dt: " << solver.dt() << "\n";
}

template<class Solver, class... Solvers>
void write_timestep(const Solver& solver, const Solvers&... solvers) noexcept {
  write_timestep(solver);
  write_timestep(solvers...);
}

template<SInd nd, class Solver>
void run_solver
(grid::Grid<nd>& grid, Solver& solver,
 const Ind maxNoTimeSteps, const Ind outputInterval) {
  initialize(grid, solver);
  write_domains(grid, solver);
  while (!solver_finished(maxNoTimeSteps, solver)) {
    write_timestep(solver);
    solver.solve();
    write_output(outputInterval, solver);
    if (solution_diverged(solver)) {
      write_domains(solver);
      TERMINATE("Solution diverged!");
    }
  }
  write_domain(solver);
}

////////////////////////////////////////////////////////////////////////////////
} // namespace fv
} // namespace solver
} // namespace hom3
////////////////////////////////////////////////////////////////////////////////
#endif
