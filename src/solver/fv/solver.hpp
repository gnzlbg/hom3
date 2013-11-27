#ifndef HOM3_SOLVERS_FV_SOLVER_HPP_
#define HOM3_SOLVERS_FV_SOLVER_HPP_
////////////////////////////////////////////////////////////////////////////////
/// Includes:
#include <limits>
#include <vector>
#include <algorithm>
#include "grid/grid.hpp"
#include "solver/fv/boundary_condition.hpp"
#include "solver/fv/container.hpp"
#include "solver/fv/tags.hpp"
#include "geometry/algorithms.hpp"
#include "quadrature/quadrature.hpp"
/// Options:
#define ENABLE_DBG_ 0
#include "misc/dbg.hpp"
////////////////////////////////////////////////////////////////////////////////
namespace hom3 {
////////////////////////////////////////////////////////////////////////////////

/// \brief Hom3 solvers
namespace solver {

/// \brief Hom3 Finite Volume Solver
namespace fv {

/// \brief Checks that \p lIdx_ and \p rIdx_ are neighbors located at opposite
/// positions of each other.
#define assert_opposite_neighbors(lIdx_, rIdx_)                         \
  using std::to_string;                                                 \
  ASSERT(is_valid((lIdx_)), "invalid lIdx!");                           \
  ASSERT(is_valid((rIdx_)), "invalid rIdx!");                           \
  ASSERT(is_valid(which_neighbor((lIdx_), (rIdx_))),                    \
         "lIdx: " + to_string((lIdx_)) + " and rIdx: "                  \
         + to_string((rIdx_)) + " aren't neighbors!");                  \
  ASSERT(is_valid(which_neighbor((rIdx_), (lIdx_))),                    \
         "lIdx: " + to_string((lIdx_)) + " and rIdx: "                  \
         + to_string((rIdx_)) + " aren't neighbors!");                  \
  ASSERT(which_neighbor(lIdx_, rIdx_)                                   \
         == grid().                                                     \
         opposite_neighbor_position(which_neighbor(rIdx_, lIdx_)),      \
         "lIdx: " + to_string((lIdx_)) + " and rIdx: "                  \
         + to_string((rIdx_)) + " aren't opposite neighbors!")

/// \brief Finite Volume solver
///
/// Implements a finite volume solver that is:
/// - first-order accurate for hyperbolic terms, and
/// - second-order accurate for parabollic terms.
///
/// Requirements on Physics component
/// - template<class T> auto physics_output(T& output_stream) const;
/// - template<class _> NumA<nvars> compute_num_flux
///   (CellIdx lIdx, CellIdx rIdx, SInd d, Num dx = opt, Num dt = opt) const;
/// - template<class _> NumA<nvars> compute_source_term(CellIdx cIdx) const;
/// - template<class _> Num compute_dt(CellIdx cIdx) const;
///
/// Optional requirements on physics component:
/// - template<class _> bool check_variables(CellIdx cIdx) const;
template<template <class> class PhysicsTT, class TimeIntegration>
struct Solver : PhysicsTT<Solver<PhysicsTT, TimeIntegration>> {
  /// \name Type traits
  ///@{
  using Physics           = PhysicsTT<Solver<PhysicsTT, TimeIntegration>>;
  using physics_type      = typename Physics::physics_type;
  static const SInd nd    = Physics::nd;
  static const SInd nvars = Physics::nvars;

  using CellContainer     = Container<nd, nvars>;
  using InitialCondition  = std::function<NumA<nvars>(const NumA<nd>)>;
  using InitialDomain     = std::function<bool(const NumA<nd>)>;

  using Boundary          = bc::Interface<nd>;
  using Boundaries        = std::vector<Boundary>;

  using Grid              = grid::CartesianHSP<nd>;
  ///@}

  /// \name Solver interface
  ///@{

  /// \brief Solver constructor
  ///
  /// Required properties are:
  /// - grid
  /// - maxNoCells
  /// - any extra properties required by the Physics class
  Solver(SolverIdx solverId, io::Properties input)
    : Physics(input)
    , solverIdx_(SolverIdx{solverId})
    , properties_(input)
    , grid_(*(io::read<Grid*>(input, "grid")))
    , cells_(io::read<Ind>(input, "maxNoCells"))
    , firstGC_(invalid<CellIdx>())
    {}
  ~Solver() {}

  /// \brief Returns the solver id
  /// \complexity O(1)
  inline SolverIdx solver_idx() const noexcept { return solverIdx_; }
  /// \brief Returns the current solution step
  /// \complexity O(1)
  inline Ind step() const noexcept { return step_; }
  /// \brief Returns the current dimensionless solution time
  /// \complexity O(1)
  inline Num time() const noexcept { return time_; }
  /// \brief Returns the min required time-step
  /// \complexity O(N)
  inline Num min_dt() const noexcept { return compute_dt(); }
  /// \brief Returns the current dimensionless solution time-step
  /// \complexity O(1)
  inline Num dt() const noexcept { return dt_; }
  /// \brief Specifies the dimensionless solution time-step to be used
  /// \complexity O(1)
  inline void dt(const Num dt__) noexcept { return force_dt(dt__); }
  /// \brief Returns the maximum dimensionless solution time allowed
  /// \warning executing solve for time >= final_time is undefined
  /// \complexity O(1)
  inline Num final_time() const noexcept { return tEnd_; }
  /// \brief Maps a local id to a global id
  /// \complexity O(1)
  const NodeIdx node_idx(const CellIdx cIdx) const noexcept {
    return cells().node_idx(cIdx);
  }

  /// \brief Returns the solver type-name (string)
  /// \warning for I\O only!
  /// For boolean operations please use solver_type() (unimplemented?)
  /// \complexity O(1)
  static inline String solver_type_name() noexcept { return "Fv"; }

  /// \brief Initializes the solver
  void initialize() {
    init_solver_variables();
    create_local_cells();
    impose_initial_condition();
    set_dt();
  }
  /// \brief Advances the solution by a single step
  void solve() noexcept {
    apply_bcs(lhs);
    set_dt();
    evolve();
    advance();
  }
  ///@}

  // Following members are public but not part of the solver interface

  /// \name Input
  ///@{
  /// \brief Sets an initial condition
  void set_initial_condition(InitialCondition initialCondition) {
    io::insert<InitialCondition>
      (properties_, "initialCondition", initialCondition);
  }

  /// \brief Appends a boundary condition
  void append_bc(Boundary bc) {
    boundaryConditions_.push_back(bc);
    grid().append_boundary(bc);
  }

  /// \brief Appends boundary conditions
  void append_bcs(Boundaries bcs) {
    for (auto&& bc : bcs) {
      append_bc(bc);
    }
  }

  ///@}

  /// \name Grid/Cell data accessors/ranges
  ///@{

  /// \brief Grid accessors
  inline const Grid& grid() const noexcept { return grid_; }
  inline       Grid& grid()       noexcept { return grid_; }

  /// \brief Cell accessors
  inline       CellContainer& cells()       noexcept { return cells_; }
  inline const CellContainer& cells() const noexcept { return cells_; }

  /// \brief Range over all variables
  inline Range<SInd> variables() const noexcept
  { return {SInd{0}, SInd{nvars}}; }

  /// \brief Conservative variables
  template<class T>
  inline NumA<nvars> Q(const CellIdx cIdx) const noexcept
  { return Q(T(), cIdx); }
  template<class T>
  inline Eigen::Block<NumM<nvars>, 1, nvars> Q(const CellIdx cIdx) noexcept
  { return Q(T(), cIdx); }
  template<class T>
  inline Num  Q(const CellIdx cIdx, const SInd varId) const noexcept
  { return Q(T(), cIdx, varId); }
  template<class T>
  inline Num& Q(const CellIdx cIdx, const SInd varId) noexcept
  { return Q(T(), cIdx, varId); }

  /// \name Tag-dispatching for conservative variables (rhs/lhs)
  ///@{
  inline  Num  Q(rhs_tag, const CellIdx cIdx, const SInd varId) const noexcept
  { return cells().rhs(cIdx, varId); }
  inline  Num& Q(rhs_tag, const CellIdx cIdx, const SInd varId)       noexcept
  { return cells().rhs(cIdx, varId); }
  inline NumA<nvars> Q(rhs_tag, const CellIdx cIdx) const noexcept
  { return cells().rhs.row(cIdx); }
  inline Eigen::Block<NumM<nvars>, 1, nvars> Q
  (rhs_tag, const CellIdx cIdx) noexcept
  { return cells().rhs.row(cIdx); }
  inline Num  Q(lhs_tag, const CellIdx cIdx, const SInd varId) const noexcept
  { return cells().lhs(cIdx, varId); }
  inline Num& Q(lhs_tag, const CellIdx cIdx, const SInd varId)       noexcept
  { return cells().lhs(cIdx, varId); }
  inline NumA<nvars> Q(lhs_tag, const CellIdx cIdx) const noexcept
  { return cells().lhs.row(cIdx); }
  inline Eigen::Block<NumM<nvars>, 1, nvars> Q
  (lhs_tag, const CellIdx cIdx) noexcept
  { return cells().lhs.row(cIdx); }
  ///@}
 private:
  /// \todo C&P code: REFACTOR: see grid/container.hpp (move to grid/range.hpp?)
  /// \brief Returns [RangeFilter] of existing node ids

  inline RangeFilter<CellIdx> valid() const noexcept
  { return {[&](const CellIdx cIdx) { return is_valid(cIdx); }}; }

  inline RangeFilter<CellIdx> ghosts() const noexcept
  { return {[&](const CellIdx cIdx) { return is_ghost_cell(cIdx); }}; }

  inline RangeFilter<CellIdx> not_ghosts() const noexcept
  { return {[&](const CellIdx cIdx) { return !is_ghost_cell(cIdx); }}; }

  inline RangeFilter<CellIdx> boundary() const noexcept
  { return {[&](const CellIdx cIdx) {
        return is_valid(cells().bc_idx(cIdx)); }}; }

  inline RangeFilter<CellIdx> boundary(const SInd bcIdx) const noexcept
  { return {[&](const CellIdx cIdx) {
        return cells().bc_idx(cIdx) == bcIdx; }}; }

  inline RangeTransformer<CellIdx, NodeIdx> cell_to_node() const noexcept
  { return {[&](const CellIdx cIdx) { return node_idx(cIdx); }}; }

  /// \todo to remove when auto works again with lambdas
  /// (see all_neighbors below: the problem is with RETURNS you can't use
  /// a lambda directly
  inline auto pos_to_neighbor(const CellIdx cIdx) const noexcept
  -> RangeTransformer<SInd, CellIdx> {
    return { [&, cIdx](const SInd nghbrPos) {
        return cells().neighbors(cIdx, nghbrPos);
      } };
  }

  /// \brief Range of all local ids (previously all_cells)
  /// \warning This includes local ids without a global id
  inline auto cell_ids() const RETURNS(cells().all_cells());

  inline auto ghost_cells() const RETURNS(cell_ids() | ghosts());
  inline auto boundary_cells() const RETURNS(cell_ids() | boundary());
  inline auto boundary_cells(const SInd bcIdx) const
  RETURNS(cell_ids() | boundary(bcIdx));

  inline auto all_neighbors(const CellIdx cIdx) const
  RETURNS(grid().neighbor_positions() | pos_to_neighbor(cIdx));

  inline auto neighbors(const CellIdx cIdx) const
  RETURNS(all_neighbors(cIdx) | valid());

 public:
  /// \brief Range of all internal cells
  /// i.e. all cells which are not ghost cells
  /// \warning internal_cells assumes that the Ghost Cells have
  /// been sorted!
  inline auto internal_cells() const -> Range<CellIdx> {
    return boost::counting_range(CellIdx{0}, is_valid(firstGC_)?
                                 firstGC_ : cells().last());
  }

  /// \brief Range of global ids
  inline auto node_ids() const
  RETURNS(cell_ids() | cell_to_node() | grid().valid());

  /// \brief Boudnary information: cell indices and positions wrt each other.
  struct BndryInfo {
    const CellIdx bndryIdx;  ///< Boundary cell index.
    const SInd    bndryPos;  ///< Boundary cell position wrt the ghost cell.
    const CellIdx ghostIdx;  ///< Ghost cell index.
    const SInd    ghostPos;  ///< Ghost cell position wrt the boundary cell.
  };

  /// \brief Returns the ids and positions of the cells involved in the boundary
  /// associated with the ghost cell \p ghostCellIdx
  ///
  /// Note: ghost cells can only have one boundary cell associated with them.
  inline auto boundary_info(const CellIdx ghostCellIdx) const noexcept
  -> BndryInfo {
    ASSERT(is_ghost_cell(ghostCellIdx), "Expecting a ghost cell!");
    SInd bndryPos = invalid<SInd>();
    CellIdx bndryIdx = invalid<CellIdx>();
    for (const auto nghbrPos : grid().neighbor_positions()) {
      if (is_valid(cells().neighbors(ghostCellIdx, nghbrPos))) {
        bndryPos = nghbrPos;
        bndryIdx = cells().neighbors(ghostCellIdx, nghbrPos);
        break;
      }
    }
    ASSERT(is_valid(bndryIdx), "boundary cell not found!");
    ASSERT(bndryPos != invalid<SInd>(), "boundary cell not found!");
    const SInd ghostPos = grid().opposite_neighbor_position(bndryPos);
    return {bndryIdx, bndryPos, ghostCellIdx, ghostPos};
  }

  ///@}


  /// \name Output
  ///@{

  /// \brief Returns the domain name
  inline String domain_name() const noexcept {
    return solver_type_name() + "_"
        + this->physics_name() + "_"
        + "Id" + to_string(solver_idx());
  }

  /// \brief Writes solver domain to VTK
  friend void write_domain(Solver& solver) noexcept {
    using std::to_string;
    String fName = solver.domain_name() + "_" + to_string(solver.step());
    solver.apply_bcs(lhs);
    write_domain(fName, solver);
  }

  /// \brief Writes solver domain to VTK
  friend void write_domain(const String fName, const Solver& solver) noexcept {
    std::cerr << "Writing domain: " << solver.domain_name() << " "
              << "| Step: " << solver.step() << " "
              << "| Time: " << solver.time() << "\n";
    io::Vtk<nd, io::format::binary> out
      (io::StreamableDomain<nd>
        (  // the following is not necessary but really cool:
          [&]() -> AnyRange<Ind> {
            return {algorithm::join(solver.cell_ids() | solver.not_ghosts(),
                                    solver.cell_ids() | solver.ghosts())
                   | boost::adaptors::transformed([](CellIdx i) { return i(); })
          }; },
          [&](const Ind cIdx) { return solver.cell_vertices(CellIdx{cIdx}); }),
        fName, io::precision::standard());

    out << io::stream("locallCellIds", 1, [](const Ind cId, const SInd) {
      return cId;
    });
    out << io::stream("gridNodeIds", 1, [&](const Ind cId, const SInd) {
      return solver.node_idx(CellIdx{cId})();
    });
    out << io::stream
        ("local_nghbrs", solver.grid().no_samelvl_neighbor_positions(),
         [&](const Ind cIdx, const SInd pos) -> Ind {
           return solver.cells().neighbors(CellIdx{cIdx}, pos)();
    });
    out << io::stream("ghost_cell", 1, [&](const Ind cIdx, const SInd) -> Ind {
      return solver.is_ghost_cell(CellIdx{cIdx}) ? 1 : invalid<Ind>();
    });
    out << io::stream("bcIdx", 1, [&](const Ind cIdx, const SInd) -> Ind {
      return solver.cells().bc_idx(CellIdx{cIdx});
    });

    out << io::stream(Solver::V::cv_names, nvars,
      [&](const Ind cIdx, const SInd d) -> Num {
        return solver.cells().lhs(CellIdx{cIdx}, d);
    });

    solver.physics()->template physics_output(out);
  }
  ///@}

 private:
  /// \name General data
  ///@{

  /// Stores the solver id
  const SolverIdx solverIdx_;
  /// Solver properties
  io::Properties properties_;
  /// Solver grid
  Grid& grid_;
  /// Cell container
  CellContainer cells_;
  /// Boundary conditions
  Boundaries boundaryConditions_;

  Boundaries& boundary_conditions() noexcept { return boundaryConditions_; }
  const Boundaries& boundary_conditions() const noexcept
  { return boundaryConditions_; }
  const Boundary& boundary_condition(const SInd bcIdx) const noexcept {
    return boundaryConditions_[bcIdx];
  }

  ///@}

  /// \name Numerical data
  ///@{
  /// Time-step
  Num dt_;
  /// Solution time
  Num time_;
  /// Max allowed solution time
  Num tEnd_;
  /// Solution step
  Ind step_;
  /// Forced time-step
  Num forced_dt_;
  /// Step at which the dt is to be forced to equal forced_dt_
  Ind forced_dt_step_;

  ///@}

  /// Offset to the first ghost cell in the container
  CellIdx firstGC_;

  /// Offset to the last ghost cell in the container
  CellIdx lastGC_;

  // Avoids insanity while working with Physics CRTP
  friend Physics;

  /// \brief Access to Physics class
  inline       Physics* physics()       noexcept
  { return static_cast<Physics*>(this); }
  inline const Physics* physics() const noexcept
  { return static_cast<const Physics*>(this); }

  /// \name Numerical functions
  ///@{

  /// \brief Computes the numerical flux
  template<class T>
  inline NumA<nvars> num_flux(const CellIdx cIdx) const noexcept {
    DBG("first term in rhs:");
    DBGV((cIdx));
    NumA<nvars> result = NumA<nvars>::Zero();
    const auto dx = cells().length(cIdx);
    for (auto d : grid().dimensions()) {
      const SInd nghbrM = d * 2;
      const SInd nghbrP = nghbrM + 1;
      const auto nghbrMId = cells().neighbors(cIdx, nghbrM);
      const auto nghbrPId = cells().neighbors(cIdx, nghbrP);
      const auto flux_m = physics()->template
                          compute_num_flux<T>(nghbrMId, cIdx, d, dx, dt());
      const auto flux_p = physics()->template
                          compute_num_flux<T>(cIdx, nghbrPId, d, dx, dt());
      result += (flux_m - flux_p);
      DBGV((cIdx)(d)(dx)(dt())(nghbrMId)(nghbrPId)(flux_m)(flux_p)(result)
           (Q<T>(nghbrMId))(Q<T>(cIdx))(Q<T>(nghbrPId)));
    }
    // return dt() / V * A * result.array();
    return dt() / dx * result.array();
  }

  template<class T>
  inline NumA<nvars> source_term(T, const CellIdx cIdx) const noexcept {
    return physics()->template compute_source_term(T(), cIdx);
  }

  /// \brief Performs an 1st-order Euler-Forward step for all cells in range \p
  /// cells
  template<class CellIdxRange>
  inline void evolve(CellIdxRange&& cells,
                     time_integration::euler_forward) noexcept {
    for (auto&& cIdx : cells) {
      Q(rhs, cIdx) = Q(lhs, cIdx) + num_flux<lhs_tag>(cIdx).transpose()
                     + source_term(lhs, cIdx).transpose();
    }
    for (auto&& cIdx : cells) {
      Q(lhs, cIdx) = Q(rhs, cIdx);
    }
  }

  /// \brief Performs a 2nd-order RK2 step for all cells in range \p cells
  template<class CellIdxRange>
  inline void evolve(CellIdxRange&& cells,
                     time_integration::runge_kutta_2) noexcept {
    for (auto&& cIdx : cells) {
      Q(rhs, cIdx) = Q(lhs, cIdx) + num_flux<lhs_tag>(cIdx).transpose();
    }

    apply_bcs(rhs);
    for (auto&& cIdx : cells) {
      Q(lhs, cIdx) = 0.5 * (Q(lhs, cIdx) + Q<rhs_tag>(cIdx)
                           + num_flux<rhs_tag>(cIdx).transpose());
    }
  }

  /// \brief Integrates the solution in time
  inline void evolve() noexcept {
    evolve(internal_cells(), TimeIntegration());
  }

  /// \brief Advances the solution to the next time-step
  inline void advance() noexcept {
    time_ += dt();
    ++step_;
  }

  /// \brief Computes the time-step
  ///
  /// Computes the time-step dt as the minimum cell time-step over all internal
  /// cells
  inline Num compute_dt() const noexcept {
    auto dt = std::numeric_limits<Num>::max();
    for (auto&& cIdx : internal_cells())  {
      dt = std::min(dt, physics()->template compute_dt<lhs_tag>(cIdx));
    }
    return dt;
  }

  /// \brief Use the time-step \p dtForce for the current solution step
  ///
  /// Note: the current solution step is the result of step()
  inline void force_dt(const Num dtForced) noexcept {
    forced_dt_step_ = step();
    forced_dt_ = dtForced;
  }

  /// \brief Returns a solution step in which instead of computing the
  /// time-step (using compute_dt), its value will be forced.
  ///
  /// Note: The time-step will be forced to whatever value is passed
  /// via force_dt before calling solve.
  /// \warning Forcing the time-step without setting it to a value
  /// using force_dt results in undefined behavior!
  inline Ind forced_dt_step() const noexcept { return forced_dt_step_; }

  /// \brief Sets the solver time-step
  ///
  /// Note: when the solution is close to final_time
  /// the time-step dt is adjusted to arrive
  /// exactly at final_time
  inline void set_dt() noexcept {
    if (time_ + dt() < tEnd_) {
      if (step() == forced_dt_step()) {
        dt_ = forced_dt_;
      } else {
        dt_ = compute_dt();
      }
    } else {
      dt_ = tEnd_ - time_;
    }
  }

  /// \brief Applies all boundary conditions
  // better: sort ghost cells by bcIdx
  // then get a [fromGhostCell,toGhostCell) range
  // apply bcIdx kernel to range
  //
  // unsolved: kernel needs to access the boudnary cell
  // this access is random access
  // it would be nice to eliminate this random access
  template<class _> void apply_bcs(_) noexcept {
    namespace csa = container::sequential::algorithm;
    auto firstBndryCell
      = csa::find_if(cells(), [&](const CellIdx i) {
          return cells().bc_idx(i) != invalid<SInd>();
    });
    // this should check the #of bndry conditions for the solver,
    // if it is indeed 0 nothing should happen, but if its not, ASSERT is fine
    ASSERT(firstBndryCell != cells().last(), "no boundary cells found!");
    auto bcIdx = cells().bc_idx(firstBndryCell);  // first bcIdx
    auto eql_bcIdx = [&](const CellIdx i) {
      return cells().bc_idx(i) == bcIdx;
    };
    auto firstGhostCell = csa::find_if(cells(), eql_bcIdx);
    auto lastGhostCell = firstGhostCell;
    // DBGV((bcIdx)(firstGhostCell));
    for (auto& boundaryCondition : boundary_conditions()) {
      ++bcIdx;
      firstGhostCell = lastGhostCell;
      lastGhostCell = csa::find_if(firstGhostCell, cells().last(), eql_bcIdx);
      // DBGV((physics()->physics_name())(bcIdx)
      // (firstGhostCell)(lastGhostCell));
      auto GCRange = Range<CellIdx>{firstGhostCell, lastGhostCell};
      boundaryCondition.apply(_(), GCRange);
    }
  }

 public:
  /// \brief Distance between cell centers of \p a and \p b
  auto cell_dx(CellIdx a, CellIdx b)
  RETURNS(geometry::algorithm::distance(cells().x_center.row(a).transpose(),
                                        cells().x_center.row(b).transpose()));

  /// \brief Computes the slope of the variable \p v at the center of cell \p
  /// cIdx in direction \p dir
  template<class _>
  auto slope(const CellIdx cIdx, const SInd v, const SInd dir) const noexcept {
    using namespace container::hierarchical;  // todo remove!
    const auto nghbrNegIdx
      = cells().neighbors(cIdx, neighbor_position(dir, neg_dir));
    const auto nghbrPosIdx
      = cells().neighbors(cIdx, neighbor_position(dir, pos_dir));

    if (is_valid(nghbrNegIdx) && is_valid(nghbrPosIdx)) {
      // both neighbors exists
      return (Q(_(), nghbrPosIdx, v) - Q(_(), nghbrNegIdx, v))
          / (2. * cells().length(cIdx));
    } else if (is_valid(nghbrNegIdx)) {  // only neg neighbor exists
      return slope<_>(nghbrNegIdx, v, dir);
    } else if (is_valid(nghbrPosIdx)) {  // only pos neighbor exists
      return slope<_>(nghbrPosIdx, v, dir);
    } else {
      // no neighbors exist in direction dir == cIdx is a ghostCell
      ASSERT(no_nghbrs(cIdx) == 1, "Ghost cells should have one neighbor!");
      ASSERT(std::begin(neighbors(cIdx)) != std::end(neighbors(cIdx)),
             "ghost cell with no neighbors?");
      ASSERT(is_ghost_cell(cIdx), "cIdx is not a ghostCell ?!?!");

      const auto bndryIdx = boundary_info(cIdx).bndryIdx;
      const auto bcIdx = cells().bc_idx(cIdx);
      return boundary_condition(bcIdx).slope(_(), bndryIdx, v, dir);
    }
  }

  /// \brief Returns the position of \p nghbrIdx w.r.t. \p cIdx and returns
  /// an invalid position if they are not neighbors.
  SInd which_neighbor(const CellIdx cIdx,
                      const CellIdx nghbrIdx) const noexcept {
    using std::to_string;
    for (auto pos : grid().neighbor_positions()) {
      if (cells().neighbors(cIdx, pos) == nghbrIdx) {
        return pos;
      }
    }
    return invalid<SInd>();
  }

  /// \brief Conservative variable \p v at surface between cells \p lIdx and \p
  /// rIdx
  template<class _> auto Q_srfc
  (_, const CellIdx lIdx, const CellIdx rIdx, const SInd v) const noexcept {
    assert_opposite_neighbors(lIdx, rIdx);
    return at_surface(lIdx, rIdx, [&, v](const CellIdx cIdx) {
        return Q(_(), cIdx, v); });
  }

  /// \brief Evaluates \p f at the surface between cells \p lIdx and \p rIdx
  template<class F> auto at_surface
  (const CellIdx lIdx, const CellIdx rIdx, F&& f) const noexcept {
    assert_opposite_neighbors(lIdx, rIdx);
    return 0.5 * (f(lIdx) + f(rIdx));
  }

  /// \brief Computes the slope of the variable \p var at the surface between
  /// cells \p lIdx and \p rIdx in direction \p slopeDir
  template<class _>
  auto surface_slope(_, const CellIdx lIdx, const CellIdx rIdx,
                     const SInd v, const SInd slopeDir) const noexcept {
    using container::hierarchical::neighbor_direction;
    const auto rIdxPosWrtLIdx = which_neighbor(lIdx, rIdx);
    ASSERT(is_valid(rIdxPosWrtLIdx), "lIdx = " + to_string(lIdx)
           + " and rIdx = " + to_string(rIdx) + " are not neighbors!");

    const auto nghbrDir = neighbor_direction(rIdxPosWrtLIdx);
    ASSERT(cells().x_center(rIdx, nghbrDir) > cells().x_center(lIdx, nghbrDir),
           "cells are in the wrong order!");

    if (nghbrDir == slopeDir) {
      // if they are neighbors in the same direction as slopeDir then use a
      // central difference:
      auto dx = cells().length(rIdx);
      ASSERT(math::approx(dx, cells().length(lIdx)), "different cell lengths!");
      return (Q(_(), rIdx, v) - Q(_(), lIdx, v)) / dx;
    } else {
      // otherwise: average the slopes
      return 0.5 * (slope<_>(rIdx, v, slopeDir) + slope<_>(lIdx, v, slopeDir));
    }
  }

  ///@}

  private:
  /// \name Debugging utilities
  ///@{

  /// \brief Checks that cell variables satisfy
  /// meaninguf physical requirements as specified
  /// by Physics::check_variables(cellId)
  void check_variables() const noexcept {
    for (auto cId : internal_cells()) {
      physics()->check_variables(cId);
    }
  }

  /// \brief Checks that all global nghbrs of local cell \p cIdx agree with
  /// those in the grid
  void check_nghbrs(const CellIdx cIdx) const noexcept {
    const auto nodeIdx = node_idx(cIdx);
    ASSERT(is_valid(nodeIdx),
           "you can only check cells that have valid global ids!");
    const auto nodeNghbrIds
      = grid().template all_samelvl_neighbors<strict>(nodeIdx);
    SInd nghbrPos = 0;
    for (auto nodeNghbrIdx : nodeNghbrIds) {
      // is there a global neighbor there belonging to the same domain?
      if (is_valid(nodeNghbrIdx)
         && grid().has_solver(nodeNghbrIdx, solver_idx())) {
        const auto cellNghbrIdx = grid().cell_idx(nodeNghbrIdx, solver_idx());
        ASSERT(cells().neighbors(cIdx, nghbrPos) == cellNghbrIdx,
               "cell: cIdx = " << cIdx << " (nodeIdx = " << node_idx(cIdx)
               <<  ") has a wrong nghbr in pos " << nghbrPos
               << "(gridId: " << cellNghbrIdx << ", local solverId: "
               << cells().neighbors(cIdx, nghbrPos));
      }
      ++nghbrPos;
    }
  }

  /// \brief Checks neighbros for all internal cells
  bool check_all_nghbrs() const noexcept {
    for (auto cIdx : internal_cells()) {
      check_nghbrs(cIdx);
    }
    return true;
  }

  /// \brief Checks the grid/solver links of cell \p cIdx
  void check_cell(const CellIdx cIdx) const noexcept {
    const auto nodeIdx = node_idx(cIdx);
    ASSERT(is_valid(nodeIdx),
           "you can only check cells that have valid global ids!");
    ASSERT(grid().cell_idx(nodeIdx, solver_idx()) == cIdx,
           "grid node has wrong cellIdx!");
  }

  /// \brief Checks the grid/solver links of all internal cells
  bool check_all_cells() const noexcept {
    for (auto cIdx : internal_cells()) {
      check_cell(cIdx);
    }
    return true;
  }

  ///@}

  /// \name Grid functions
  ///@{

  /// \brief Number of cell neighbors (including ghost-cells)
  SInd no_nghbrs(const CellIdx cIdx) const noexcept {
    SInd noNghbrs = 0;
    for (const auto nghbrPos : grid().neighbor_positions()) {
      if (is_valid(cells().neighbors(cIdx, nghbrPos))) {
        ++noNghbrs;
      }
    }
    return noNghbrs;
  }

  /// \brief Lazy range of all cell vertices
  auto cell_vertices(const CellIdx cIdx) const noexcept
  -> typename Grid::CellVertices {
    const auto x_c = cells().x_center.row(cIdx);

    const auto length
        = !is_ghost_cell(cIdx)
        // cell is a global cell -> get its length
        ? grid().cell_length(node_idx(cIdx))
        // cell is a ghost cell -> get its boundary cell's length
        : [&]() {
            const auto bndryCellId = boundary_info(cIdx).bndryIdx;
            return grid().cell_length(node_idx(bndryCellId));
    }();

    return {cIdx(), grid().cell_vertices_coords(length, x_c)};
  }
  ///@}

  /// \name Ghost-cell related
  ///@{

  inline bool is_ghost_cell(const CellIdx cIdx) const noexcept
  {
    ASSERT([&]() {
        if(!is_valid(node_idx(cIdx))) {
          ASSERT(is_valid(cells().bc_idx(cIdx)), "Error: ghost cell does not have bcIdx!");
        }
        return true; }(), "ghost cell check!");
    return !is_valid(node_idx(cIdx)) ; }

  /// \brief Computes the cell-center coordinates of the ghost cell \p
  /// ghostCellIdx
  NumA<nd> ghost_cell_coordinates(const CellIdx ghostCellIdx) const noexcept  {
    ASSERT(no_nghbrs(ghostCellIdx) == 1, "Invalid ghost cell!");

    /// Find local boundary cell id and the ghost cell position w.r.t. the
    /// boundary cell
    auto const bndryInfo = boundary_info(ghostCellIdx);
    auto const bndryCellIdx = bndryInfo.bndryIdx;
    auto pos_wrt_ghostCell = bndryInfo.bndryPos;
    ASSERT(is_valid(pos_wrt_ghostCell), "cells are not neighbors!");
    ASSERT(is_valid(bndryCellIdx), "ghost cell has no associated bndry cell!");

    /// \todo: opposite nghbr position should be a non-member function!
    const auto pos_wrt_boundary
        = grid().opposite_neighbor_position(pos_wrt_ghostCell);
    ASSERT(is_valid(pos_wrt_boundary), "cells are not neighbors!");

    auto length = cells().length(bndryCellIdx);
    auto x_gc = cells().x_center.row(bndryCellIdx).transpose() + length
                * grid().nghbr_rel_pos(pos_wrt_boundary).template cast<Num>();
    return x_gc;
  }

  /// \brief Sorts ghost cells by boundary condition id and returns the id of
  /// the first ghost cell
  CellIdx sort_gc() noexcept {
    /// Sort ghost cells by bcIdx:
    auto valid_bcIdx = [&](CellIdx i) { return is_valid(cells().bc_idx(i)); };
    const auto firstGhostCell
        = container::sequential::algorithm::find_if(cells(), valid_bcIdx);

    auto sort_bcIds = [&](typename CellContainer::value_type a,
                          typename CellContainer::value_type b) {
      return a.bc_idx() < b.bc_idx();
    };
    std::sort(cells().begin() + firstGhostCell(), cells().end(), sort_bcIds);

    /// Correct nghbrs in boundary cells:
    for (auto gcIdx : ghost_cells()) {
      for (auto nghbrPos : grid().neighbor_positions()) {
        auto nghbrIdx = cells().neighbors(gcIdx, nghbrPos);
        if (is_valid(nghbrIdx)) {
          auto oppositeNghbrPos = grid().opposite_neighbor_position(nghbrPos);
          cells().neighbors(nghbrIdx, oppositeNghbrPos) = gcIdx;
        }
      }
    }
    return firstGhostCell;
  }

  /// Creates a ghost-cell for \p localBCellId located in the position of the
  /// neighbor \p nghbrPos (w.r.t the boundary cell) and sets the ghost cell
  /// coordinates.
  ///
  /// \returns the local ghostCellIdx
  ///
  /// Precondition: cell at \p localBCellId has no neighbor in \nghbrPos
  /// Postcondition:
  ///   - ghost cell that has a single neighbor
  ///   - and a boundary condition index
  CellIdx create_ghost_cell(const CellIdx bndryCellIdx,
                            const SInd nghbrPos, const SInd bcIdx) noexcept {
    ASSERT(!is_valid(cells().neighbors(bndryCellIdx, nghbrPos)),
           "There is already a cell there!");

    /// 0) create a ghost cell without a globalBCellId
    const auto ghostCellIdx = cells().push_cell();
    cells().node_idx(ghostCellIdx) = invalid<NodeIdx>();

    /// 1) Set ghost cell boundary condition:
    cells().bc_idx(ghostCellIdx) = bcIdx;

    /// 2) Set neighbor relationships:
    const auto oppositeNghbrPos = grid().opposite_neighbor_position(nghbrPos);
    cells().neighbors(bndryCellIdx, nghbrPos) = ghostCellIdx;
    cells().neighbors(ghostCellIdx, oppositeNghbrPos) = bndryCellIdx;

    /// 3) Set ghost cell position
    cells().x_center.row(ghostCellIdx) = ghost_cell_coordinates(ghostCellIdx);
    cells().length(ghostCellIdx) = cells().length(bndryCellIdx);

    /// \todo Test post-conditions
    return ghostCellIdx;
  }

  ///@}

  /// \name Initialization
  ///@{
  void init_solver_variables() noexcept {
    dt_ = 0;
    time_ = 0;
    tEnd_ = io::read<Num>(properties_, "timeEnd");
    step_ = 0;
    forced_dt_ = 0;
    forced_dt_step_ = invalid<Ind>();
  }

  void create_local_cells() noexcept {
    // 1) create a local cell for each grid leaf cell
    //    set the nodeIdx in the local cells
    //    set the cIdx in the grid cells
    //    set the cell coordinates
    {
      auto initialDomain
          = io::read<InitialDomain>(properties_, "initialDomain");
      for (auto nIdx : grid().leaf_nodes()) {
        auto xc = grid().cell_coordinates(nIdx);
        if (!initialDomain(xc)) { continue; }
        auto cIdx = cells().push_cell();
        cells().node_idx(cIdx) = nIdx;
        grid().cell_idx(nIdx, solver_idx()) = cIdx;
        cells().x_center.row(cIdx) = xc;
        cells().length(cIdx) = grid().cell_length(nIdx);
        cells().bc_idx(cIdx) = invalid<SInd>();
      }
    }

    ASSERT(check_all_cells(), "solver cells / grid node links are wrong!");

    // set local nghbr ids: (localCellId,nIdx) -> globalNghbrIds ->
    // localNghbrIds
    for (auto cIdx : cell_ids()) {
      auto nodeNghbrIds =  grid().template
                           all_samelvl_neighbors<strict>(node_idx(cIdx));
      for (auto nghbrPos : grid().neighbor_positions()) {
        auto nghbrIdx = nodeNghbrIds(nghbrPos);
        cells().neighbors(cIdx, nghbrPos)
            = is_valid(nghbrIdx) ? grid().cell_idx(nghbrIdx, solver_idx())
                                 : invalid<CellIdx>();
      }
    }

    firstGC_ = CellIdx{cells().size()};

    create_ghost_cells();

    /// sort ghost cells by boundary id
    sort_gc();

    /// set distances to nghbrs
    {
      for (auto cIdx : cell_ids()) {
        for (const auto nghbrPos : grid().neighbor_positions()) {
          auto nghbrIdx = cells().neighbors(cIdx, nghbrPos);
          if (!is_valid(nghbrIdx)) { continue; }  // \todo necessary?
          DBGV((cIdx)(nghbrIdx)(nghbrPos));
          cells().distances(cIdx, nghbrPos) = cell_dx(cIdx, nghbrIdx);
          DBGV((cells().distances(cIdx, nghbrPos)));
        }
      }
    }
    ASSERT(check_all_nghbrs(), "internal cell nghbrIds don't agree with grid!");
  }

  /// Create Ghost Cells:
  ///
  /// \warning this only works for cutoff right now
  void create_ghost_cells() noexcept {
    memory::stack::arena<SInd, 6 + 2> stackMemory;  // 6nghbrs + 2 alignment
    auto noLeafCells = cells().size();
    SInd ghostCellBoundaryIdx = 0;
    for (auto boundary : boundary_conditions()) {
      for (auto bndryNodeIdx : node_ids() | grid().cut_by_boundary(boundary)) {
        const auto bndryCellIdx
          = grid().cell_idx(bndryNodeIdx, solver_idx());
        ASSERT(is_valid(bndryCellIdx), "invalid bndryCellIdx!");
        ASSERT(node_idx(bndryCellIdx) == bndryNodeIdx,
               "solver and grid are not synchronized");
        ASSERT(is_valid(node_idx(bndryCellIdx)),
               "the global id has to be valid!");

        // find missing neighbor positions
        auto missingNghbrPositions
            = memory::stack::make<std::vector>(stackMemory);
        for (const auto nghbrPos : grid().neighbor_positions()) {
          if (cells().neighbors(bndryCellIdx, nghbrPos) == invalid<CellIdx>()) {
            missingNghbrPositions.emplace_back(nghbrPos);
          }
        }

        auto neg_distance = [&](const SInd pos) {
          return boundary.signed_distance
          (grid().neighbor_coordinates(bndryNodeIdx, pos)) < 0.;
        };

        using boost::adaptors::filtered;
        for (auto nghbrPos : missingNghbrPositions | filtered(neg_distance)) {
          create_ghost_cell(bndryCellIdx, nghbrPos, ghostCellBoundaryIdx);
        }
      }
      ++ghostCellBoundaryIdx;
    }

    auto newTotalNoCells = cells().size();
    std::cerr << "fv container | #of leafs: " << noLeafCells
              << " | #of ghosts: " << newTotalNoCells - noLeafCells
              << " | #of cells: " << newTotalNoCells << "\n";
    ASSERT(newTotalNoCells - noLeafCells > 0, "#of ghost cells is 0!");
  }

  /// \brief Imposes the initial condition on the lhs
  void impose_initial_condition() noexcept {
    auto initialCondition
      = io::read<InitialCondition>(properties_, "initialCondition");
    for (auto cIdx : internal_cells()) {
      const auto length = grid().cell_length(node_idx(cIdx));
      const auto x_c = cells().x_center.row(cIdx).transpose();
      const auto average
        = quadrature::integrate<nd>(initialCondition, x_c, length)
          / geometry::cell::cartesian::volume<nd>(length);

      cells().lhs.row(cIdx) = average;
    }
  }
  ///@}
};

#define HOM3_FV_PHYSICS_SOLVER_() \
  template<SInd nd, class NumFlux, class TimeIntegration>               \
  struct SolverBuilder {                                                \
    template<class S> using physics = Physics<nd, NumFlux, S>;          \
    using type = solver::fv::Solver<physics, TimeIntegration>;          \
  };                                                                    \
  template<SInd nd, class NumFlux, class TimeIntegration>               \
  using Solver = typename SolverBuilder<nd, NumFlux, TimeIntegration>::type

}  // namespace fv

////////////////////////////////////////////////////////////////////////////////
}  // namespace solver
}  // namespace hom3
////////////////////////////////////////////////////////////////////////////////
#undef ENABLE_DBG_
////////////////////////////////////////////////////////////////////////////////
#endif
