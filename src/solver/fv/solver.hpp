#ifndef HOM3_SOLVERS_FV_SOLVER_HPP_
#define HOM3_SOLVERS_FV_SOLVER_HPP_
////////////////////////////////////////////////////////////////////////////////
/// Includes:
#include "../../grid/grid.hpp"
#include "../../geometry/algorithms.hpp"
#include "container.hpp"
#include "../../quadrature/quadrature.hpp"
/// Options:
#define ENABLE_DBG_ 0
#include "../../misc/dbg.hpp"
////////////////////////////////////////////////////////////////////////////////
namespace hom3 {
////////////////////////////////////////////////////////////////////////////////

namespace solver { namespace fv {

/// \name Tags for rhs and lhs variables
///@{
struct rhs_tag {};
struct lhs_tag {};
///@}

/// note an option is a template template parameter of the form:
/// - template <SInd,class> class
/// numerical methods take a bunch of options:
/// Method<option...> -> Method<template <SInd,class> class...>
/// - The first option is special for now (it has to be the physics class),
/// but it walks like a duck and quacks like a duck -> it's just an option
///
/// \todo concept checking (e.g. Boost.ConceptChecking) improves this mess
// a whole lot => Try C++14::tr::concepts_lite as soon as possible
// template<SInd nd, // #of spatial dimensions
//          template <SInd,class> class Physics, // Physics class
//          template <SInd,class> class... ConfigOptions // configuration options
//          >

/// \brief Finite Volume solver
///
/// \todo Move to the lId (localId) and gId (globalId) convention for ids
/// \todo Remove all localIds localCellIds globalIds globalCellIds cId...
/// \todo cells().globalIds rename to globalId (container variables should be
/// singular if they have only one element per cell)
///
///
/// Requirements on Physics component
/// - template<class T> auto physics_output(T& output_stream) const;
/// - template<class T> NumA<nvars> compute_num_flux(Ind cell_i, Ind cell_i+1,
///                                                  SInd d, Num dx = opt,
///                                                  Num dt = opt) const;
template<SInd nd_, // pretty ugly, fix in the future! (it's getting better)
         template <SInd,class> class PhysicsTT,//, // -> class Physics
         class NumFlux_//         template <SInd,class> class NumFluxTT
         >
struct Solver : PhysicsTT<nd_,Solver<nd_,PhysicsTT,NumFlux_/*,NumFluxTT*/>> { //,
  //                NumFluxTT<nd_,Solver<nd_,PhysicsTT,NumFluxTT>>  {

  /// \name Type traits
  ///@{
  static const SInd nd    = nd_;
  using NumFlux           = NumFlux_;
  using Physics           = PhysicsTT<nd,Solver<nd,PhysicsTT,NumFlux/*,NumFluxTT*/>>;
  using physics_type      = typename Physics::type;
  static const SInd nvars = Physics::nvars;


  //  using NumFlux           = NumFluxTT<nd,Solver<nd,PhysicsTT/*,NumFluxTT*/>>;

  using CellContainer     = Container<nd,nvars>;
  using BoundaryCondition = std::function<void(Ind)>;
  using InitialCondition  = std::function<NumA<nvars>(const NumA<nd>)>;
  using InitialDomain     = std::function<bool(const NumA<nd>)>;

  using Grid              = grid::CartesianHSP<nd>;

  using rhs               = rhs_tag;
  using lhs               = lhs_tag;
  ///@}

  /// \name Solver interface
  ///@{

  /// \brief Solver constructor
  ///
  /// Required properties are:
  /// - grid
  /// - maxNoCells
  /// - any extra properties required by the Physics class
  Solver(SolverIdx solverId, io::Properties input) :
      Physics(input),
      solverIdx_(SolverIdx{solverId}),
      properties_(input),
      grid_(*io::read<Grid*>(input,"grid")),
      cells_(io::read<Ind>(input,"maxNoCells"))
      {}

  /// \brief Returns the solver id
  /// \complexity O(1)
  inline SolverIdx solver_idx() const { return solverIdx_; }
  /// \brief Returns the current solution step
  /// \complexity O(1)
  inline Ind step() const { return step_; }
  /// \brief Returns the current solution time
  /// \complexity O(1)
  inline Num time() const { return time_; }
  /// \brief Returns the min required time-step
  /// \complexity O(N)
  inline Num min_dt() const { return compute_dt(); }
  /// \brief Returns the current solution time-step
  /// \complexity O(1)
  inline Num dt() const { return dt_; }
  /// \brief Specifies the dt to be used
  /// \complexity O(1)
  inline void dt(const Num dt__) { return force_dt(dt__); }
  /// \brief Returns the maximum solution time allowed
  /// \warning executing solve for time >= final_time is undefined
  /// \complexity O(1)
  inline Num final_time() const { return tEnd_; }
  /// \brief Maps a local id to a global id
  /// \complexity O(1)
  const NodeIdx node_idx(const CellIdx cIdx) const {
    return cells().node_idx(cIdx);
  }

  /// \brief Returns the solver type-name (string)
  /// \warning for I\O only! For boolean operations please use solver_type() (unimplemented?)
  /// \complexity O(1)
  static std::string solver_type_name() { return "Fv"; }

  /// \brief Initializes the solver
  void initialize() {
    init_solver_variables();
    create_local_cells();
    set_initial_condition();
  }
  /// \brief Advances the solution by a single step
  void solve() {
    apply_bcs();
    set_dt();
    evolve();
    advance();
  }
  ///@}

  // Following members are public but not part of the solver interface

  /// \name Grid/Cell data accessors/ranges
  ///@{

  /// \brief Grid accessors
  inline const Grid& grid() const { return grid_; }
  inline       Grid& grid()       { return grid_; }

  /// \brief Cell accessors
  inline       CellContainer& cells()       { return cells_; }
  inline const CellContainer& cells() const { return cells_; }

  /// \brief Conservative variables
  template<class T> inline NumA<nvars> Q(const CellIdx lId) const {
    return vars_(T(),lId);
  }
  template<class T> inline Eigen::Block<NumM<nvars>,1,nvars> Q(const CellIdx lId) {
    return vars_(T(),lId);
  }
  template<class T> inline Num  Q(const CellIdx lId, const SInd varId) const {
    return vars_(T(), lId, varId);
  }
  template<class T> inline Num& Q(const CellIdx lId, const SInd varId)       {
    return vars_(T(), lId, varId);
  }

private:

  /// \todo C&P code: REFACTOR: see grid/container.hpp (move to grid/range.hpp?)
  /// \brief Returns [RangeFilter] of existing node ids

  inline RangeFilter<CellIdx> valid() const {
    return {[&](const CellIdx cIdx){ return is_valid(cIdx); }};
  }

  inline RangeFilter<CellIdx> ghosts() const {
    return {[&](const CellIdx cIdx){ return is_ghost_cell(cIdx); }};
  }

  inline RangeFilter<CellIdx> not_ghosts() const {
    return {[&](const CellIdx cIdx){ return !is_ghost_cell(cIdx); }};
  }

  inline RangeFilter<CellIdx> boundary() const {
    return {[&](const CellIdx cIdx){ return is_valid(cells().bc_idx(cIdx)); }};
  }

  inline RangeFilter<CellIdx> boundary(const SInd bcId) const {
    return {[&](const CellIdx cIdx){ return cells().bc_idx(cIdx) == bcId; }};
  }

  inline RangeTransformer<CellIdx,NodeIdx> cell_to_node() const {
    return {[&](const CellIdx cIdx){ return node_idx(cIdx); }};
  }

  /// \brief Range of all local ids (previously all_cells)
  /// \warning This includes local ids without a global id
  inline auto cell_ids() const -> Range<CellIdx> {
    return cells().all_cells();
  }

  inline auto ghost_cells() const -> FRange<CellIdx> {
      return cell_ids() | ghosts();
  }

  inline auto boundary_cells() const -> FRange<CellIdx> {
      return cell_ids() | boundary();
  }

  inline auto boundary_cells(const SInd bcId) const -> FRange<CellIdx> {
      return cell_ids() | boundary(bcId);
  }

public:

  /// \brief Range of all internal cells
  /// i.e. all cells which are not ghost cells
  inline auto internal_cells() const -> FRange<CellIdx> {
      return cell_ids() | not_ghosts();
  }

  /// \brief Range of global ids
  inline auto node_ids() const -> decltype(cell_ids() | cell_to_node() | grid().nodes().valid()) {
    return cell_ids() | cell_to_node() | grid().nodes().valid();
  }

  /// \brief Returns the id of the boundary cell associated with the ghost cell
  /// \p ghostCellId as well as its position w.r.t. the ghost cell
  ///
  /// Note: use std/boost::get<pos>(returned_object) to extract the information
  /// \returns at pos 0: boundary cell id (Ind)
  /// \returns at pos 1: boundary cell position (SInd)
  ///
  /// \warning You shall not rely on whatever specific type this function might
  /// return at a given moment in time!
  ///
  /// Note: ghost cells can only have one boundary cell associated with them.
  inline auto boundary_cell_id(const CellIdx ghostCellId) const -> std::tuple<CellIdx,SInd> {
    for(const auto nghbrPos : grid().nodes().nghbr_pos()) {
      if(is_valid(cells().neighbors(ghostCellId,nghbrPos))) {
        return std::make_tuple(cells().neighbors(ghostCellId,nghbrPos),nghbrPos);
      }
    }
    /// this actually doesnt check for multiple nghbrs which can happen...
    ASSERT(false,"ghostcell are required to have only one neighbor!");
    return std::make_tuple(invalid<CellIdx>(),invalid<SInd>());
  }

  ///@}


  /// \name Output
  ///@{

  /// \brief Writes solver domain to VTK
  friend void write_domain(Solver* solver) {
    using std::to_string;
    std::string fName = solver->solver_type_name() + "_"
                        + solver->physics_name() + "_"
                        + "Id" + to_string(solver->solver_idx())
                        + "_" + to_string(solver->step());
    solver->apply_bcs();
    write_domain(fName,solver);
  }

  /// \brief Writes solver domain to VTK
  friend void write_domain(const std::string fName, const Solver* solver) {

    io::Vtk<nd,io::format::binary> out(
        io::StreamableDomain<nd>( // the following is not necessary but really cool:
            [&]() -> AnyRange<Ind> {
              return {algorithm::join(solver->cell_ids() | solver->not_ghosts(),
                                      solver->cell_ids() | solver->ghosts())
                    | boost::adaptors::transformed([](CellIdx i){ return i();})
                    }; },
            [&](const Ind lId){ return solver->cell_vertices(CellIdx{lId});}),
        fName, io::precision::standard());

    out << io::stream("locallCellIds",1,[](const Ind cId, const SInd){ return cId; });
    out << io::stream("local_nghbrs",solver->grid().nodes().no_samelvl_nghbr_pos(),
                      [&](const Ind lId, const SInd pos) -> Ind {
                        return solver->cells().neighbors(CellIdx{lId},pos)();
                      });
    out << io::stream("ghost_cell",1,[&](const Ind lId, const SInd) -> Ind {
        return solver->is_ghost_cell(CellIdx{lId}) ? 1 : invalid<Ind>();
      });
    out << io::stream("bcId",1,[&](const Ind lId, const SInd) -> Ind {
        return solver->cells().bc_idx(CellIdx{lId});
      });

    out << io::stream(Solver::V::cv_names,nvars,[&](const Ind cId, const SInd d) -> Num {
        return solver->cells().lhs(CellIdx{cId}, d);
      });

    solver->physics()->template physics_output(out);
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

  ///@}

  /// \name Numerical data
  ///@{

  /// \brief Returns the cfl number
  Num cfl() const { return cfl_; }
  /// CFL number
  Num cfl_;
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
  friend PhysicsTT<nd,Solver<nd,PhysicsTT,NumFlux/*,NumFluxTT*/>>;

  /// \brief Access to Physics class
  // Note: think of it as the inverse of a CRTP
  inline       Physics* physics()       { return static_cast<Physics*>(this); }
  inline const Physics* physics() const { return static_cast<const Physics*>(this); }

  /// \name Tag-dispatching for conservative variables (rhs/lhs)
  ///@{
  inline  Num  vars_(rhs, const CellIdx lId, const SInd varId) const {
    return cells().rhs(lId,varId);
  }
  inline  Num& vars_(rhs, const CellIdx lId, const SInd varId)       {
    return cells().rhs(lId,varId);
  }
  inline NumA<nvars> vars_(rhs, const CellIdx lId) const {
    return cells().rhs.row(lId);
  }
  inline Eigen::Block<NumM<nvars>,1,nvars> vars_(rhs,const CellIdx lId) {
    return cells().rhs.row(lId);
  }
  inline Num  vars_(lhs, const CellIdx lId, const SInd varId) const {
    return cells().lhs(lId,varId);
  }
  inline Num& vars_(lhs, const CellIdx lId, const SInd varId)       {
    return cells().lhs(lId,varId);
  }
  inline NumA<nvars> vars_(lhs, const CellIdx lId) const {
    return cells().lhs.row(lId);
  }
  inline Eigen::Block<NumM<nvars>,1,nvars> vars_(lhs, const CellIdx lId) {
    return cells().lhs.row(lId);
  }
  ///@}

  /// \name Numerical functions
  ///@{

  /// \brief Computes the numerical flux
  template<class T> inline NumA<nvars> num_flux(const CellIdx cIdx) const {
    DBG("first term in rhs:");
    DBGV((cIdx));
    NumA<nvars> result = NumA<nvars>::Zero();
    const auto dx = cells().length(cIdx);
    for(SInd d = 0; d < nd; ++d) {
      const SInd nghbrM = d*2;
      const SInd nghbrP = nghbrM + 1;
      const auto nghbrMId = cells().neighbors(cIdx,nghbrM);
      const auto nghbrPId = cells().neighbors(cIdx,nghbrP);
      const auto flux_m = physics()->template compute_num_flux<T,NumFlux>(nghbrMId,cIdx,d,dx,dt());
      const auto flux_p = physics()->template compute_num_flux<T,NumFlux>(cIdx,nghbrPId,d,dx,dt());
      result += dt() / dx * (flux_m - flux_p);
      DBGV((cIdx)(d)(dx)(dt())(nghbrMId)(nghbrPId)(flux_m)(flux_p)(result)
           (Q<T>(nghbrMId))(Q<T>(cIdx))(Q<T>(nghbrPId)));
    }
    return result;
  }

  /// \brief Performs an Euler-Forward step for all cells in range \p cells
  template<class CellIdxRange>
  inline void ef(CellIdxRange&& cells) {
    apply_bcs();
    for(auto&& cIdx : cells) {
      Q<rhs>(cIdx) = Q<lhs>(cIdx) + num_flux<lhs>(cIdx).transpose();
    }
    for(auto&& cIdx : cells) {
      Q<lhs>(cIdx) = Q<rhs>(cIdx);
    }
  }

  /// \brief Performs a RK2 step for all cells in range \p cells
  ///
  /// \todo apply_bcs should be parametrized to with rhs/lhs to decide when what
  /// should be applied
  template<class CellIdxRange>
  inline void rk2(CellIdxRange&& cells) {
    apply_bcs();
    for(auto&& cIdx : cells) {
      Q<rhs>(cIdx) = Q<lhs>(cIdx) + num_flux<lhs>(cIdx).transpose();
    }

    apply_bcs();
    for(auto&& cIdx : cells) {
      Q<lhs>(cIdx) = 0.5 * (Q<lhs>(cIdx) + Q<rhs>(cIdx)
                           + num_flux<rhs>(cIdx).transpose());
    }
  }

  /// \brief Integrates cells in time
  // uses a hard-coded time integration method
  // \todo parametrize on time-integration method (using TP?)
  inline void evolve() {
    ef(internal_cells());
    //rk2(internal_cells());
    apply_bcs();
  }

  /// \brief Advances the solution to the next time-step
  inline void advance() {
    time_ += dt();
    ++step_;
  }

  /// \brief Computes the time-step
  ///
  /// Computes the time-step dt as the minimum cell time-step over all internal
  /// cells
  inline Num compute_dt() const {
    auto dt = std::numeric_limits<Num>::max();
    for(auto&& cIdx : internal_cells())  {
      dt = std::min(dt,physics()->template compute_dt<lhs>(cIdx));
    }
    return dt;
  }

  inline void force_dt(const Num dt__) {
    forced_dt_step_ = step();
    forced_dt_ = dt__;
  }

  inline Ind forced_dt_step() const { return forced_dt_step_; }

  /// \brief Sets the solver time-step
  ///
  /// Note: when the solution is close to final_time
  /// the time-step dt is adjusted to arrive
  /// exactly at final_time
  inline void set_dt() {
    if (time_ + dt() < tEnd_) {
      if(step() == forced_dt_step()) {
        dt_ = forced_dt_;
      } else {
        dt_ = compute_dt();
      }
    } else {
      dt_ = tEnd_ - time_;
    }
  }

  /// \brief Applies all boundary conditions
  // better: sort ghost cells by bcId
  // then get a [fromGhostCell,toGhostCell) range
  // apply bcId kernel to range
  //
  // unsolved: kernel needs to access the boudnary cell
  // this access is random access
  // it would be nice to eliminate this random access
  void apply_bcs() {
    auto firstBndryCell
        = container::sequential::algorithm::find_if
        (cells(),[&](CellIdx i){return cells().bc_idx(i) != invalid<SInd>(); });
    // this should check the #of bndry conditions for the solver,
    // if it is indeed 0 nothing should happen, but if its not, ASSERT is fine
    ASSERT(firstBndryCell != cells().last(), "no boundary cells found!");
    SInd bcId = cells().bc_idx(firstBndryCell); // first bcId
    auto eql_bcId = [&](CellIdx i) { return cells().bc_idx(i) == bcId; };
    auto firstGhostCell = container::sequential::algorithm::find_if(cells(),eql_bcId);
    auto lastGhostCell = firstGhostCell;
    //DBGV((bcId)(firstGhostCell));
    for(auto& boundary : grid().boundaries(solver_idx())) {
      for(auto&& condition : boundary.conditions()) {
        ++bcId;
        firstGhostCell = lastGhostCell;
        lastGhostCell = container::sequential::algorithm::find_if(firstGhostCell,cells().last(),eql_bcId);
        //DBGV((physics()->physics_name())(bcId)(firstGhostCell)(lastGhostCell));
        auto GCRange = Range<CellIdx>{firstGhostCell,lastGhostCell};
        (*condition)(GCRange);
      }
    }
  }

  ///@}

  /// \name Debugging utilities
  ///@{

  /// \brief Checks that cell variables satisfy
  /// meaninguf physical requirements as specified
  /// by Physics::check_variables(cellId)
  void check_variables() {
    for(auto cId : internal_cells()) {
      physics()->check_variables(cId);
    }
  }

  /// \brief Checks that all global nghbrs of local cell \p cIdx agree with
  /// those in the grid
  void check_nghbrs(const CellIdx cIdx) const {
    const auto nodeIdx = node_idx(cIdx);
    ASSERT(is_valid(nodeIdx), "you can only check cells that have valid global ids!");
    const auto nodeNghbrIds = grid().nodes().template all_samelvl_nghbrs<strict>(nodeIdx);
    SInd nghbrPos = 0;
    for(auto nodeNghbrIdx : nodeNghbrIds) {
      if(is_valid(nodeNghbrIdx) // is there a global neighbor there belonging to the same domain?
         && grid().nodes().has_solver(nodeNghbrIdx,solver_idx())) {
        const auto cellNghbrIdx = grid().nodes().cell_id(nodeNghbrIdx,solver_idx());
        ASSERT(cells().neighbors(cIdx,nghbrPos) == cellNghbrIdx,
               "cell: cIdx = " << cIdx << " (nodeIdx = " << node_idx(cIdx)
               <<  ") has a wrong nghbr in pos " << nghbrPos
               << "(gridId: " << cellNghbrIdx << ", local solverId: "
               << cells().neighbors(cIdx,nghbrPos));
      }
      ++nghbrPos;
    }
  }

  /// \brief Checks neighbros for all internal cells
  bool check_all_nghbrs() const {
    for(auto cIdx : internal_cells()) {
      check_nghbrs(cIdx);
    }
    return true;
  }

  /// \brief Checks the grid/solver links of cell \p cIdx
  void check_cell(const CellIdx cIdx) const {
    const auto nodeIdx = node_idx(cIdx);
    ASSERT(is_valid(nodeIdx), "you can only check cells that have valid global ids!");
    ASSERT(grid().nodes().cell_id(nodeIdx,solver_idx()) == cIdx, "grid node has wrong cellIdx!");
  }

  /// \brief Checks the grid/solver links of all internal cells
  bool check_all_cells() const {
    for(auto cIdx : internal_cells()) {
      check_cell(cIdx);
    }
    return true;
  }

  //@}

  /// \name Grid functions
  ///@{

  /// \brief Number of cell neighbors (including ghost-cells)
  SInd no_nghbrs(const CellIdx cIdx) const {
    SInd noNghbrs = 0;
    for(const auto nghbrPos : grid().nodes().nghbr_pos()) {
      if(is_valid(cells().neighbors(cIdx,nghbrPos))) {
        ++noNghbrs;
      }
    }
    return noNghbrs;
  }

  /// \brief Lazy range of all cell vertices
  auto cell_vertices(const CellIdx cIdx) const -> typename Grid::CellVertices {
    const auto x_c = cells().x_center.row(cIdx);

    const auto length
        = !is_ghost_cell(cIdx)
        ? grid().cell_length(node_idx(cIdx)) // cell is a global cell -> get its length
        : [&](){ // cell is a ghost cell -> get its boundary cell's length
      auto bndryCellId = invalid<CellIdx>();
      std::tie(bndryCellId,std::ignore) = boundary_cell_id(cIdx);
      return grid().cell_length(node_idx(bndryCellId));
    }();

    return { cIdx(), grid().cell_vertices_coords(length, x_c) };
  }
  ///@}

  /// \name Ghost-cell related
  ///@{

  inline bool is_ghost_cell(const CellIdx cIdx) const { return !is_valid(node_idx(cIdx)); }

  /// \brief Computes the cell-center coordinates of the ghost cell \p
  /// ghostCellId
  NumA<nd> ghost_cell_coordinates(const CellIdx ghostCellIdx)  {
    ASSERT(no_nghbrs(ghostCellIdx) == 1, "Invalid ghost cell!");

    /// Find local boundary cell id and the ghost cell position w.r.t. the
    /// boundary cell
    auto bndryCellIdx = invalid<CellIdx>();
    auto position_wrt_ghostCell = invalid<SInd>();
    std::tie(bndryCellIdx,position_wrt_ghostCell) = boundary_cell_id(ghostCellIdx);
    ASSERT(is_valid(position_wrt_ghostCell), "ERROR!");
    ASSERT(is_valid(bndryCellIdx),"ERROR!");

    /// \todo: opposite nghbr position should be a non-member function!
    const auto position_wrt_boundary = grid().nodes().opposite_nghbr_position(position_wrt_ghostCell);
    ASSERT(is_valid(position_wrt_boundary),"ERROR!");

    auto length = cells().length(bndryCellIdx);
    auto x_gc = cells().x_center.row(bndryCellIdx).transpose()
                + length * grid().nghbr_rel_pos(position_wrt_boundary).template cast<Num>();
    return x_gc;
  }

  /// \brief Sorts ghost cells by boundary condition id and returns the id of
  /// the first ghost cell
  CellIdx sort_gc() {

    /// Sort ghost cells by bcId:
    auto valid_bcIdx = [&](CellIdx i) { return is_valid(cells().bc_idx(i)); };
    const auto firstGhostCell
        = container::sequential::algorithm::find_if(cells(),valid_bcIdx);

    auto sort_bcIds = [&](typename CellContainer::value_type a,
                          typename CellContainer::value_type b) {
      return a.bc_idx() < b.bc_idx();
    };
    std::sort(cells().begin() + firstGhostCell(), cells().end(),sort_bcIds);

    /// Correct nghbrs in boundary cells:
    for(auto gcIdx : ghost_cells()) {
      for(auto nghbrPos : grid().nodes().nghbr_pos()) {
        auto nghbrIdx = cells().neighbors(gcIdx,nghbrPos);
        if(is_valid(nghbrIdx)) {
          auto oppositeNghbrPos = grid().nodes().opposite_nghbr_position(nghbrPos);
          cells().neighbors(nghbrIdx,oppositeNghbrPos) = gcIdx;
        }
      }
    }
    return firstGhostCell;
  }

  /// Creates a ghost-cell for \p localBCellId located in the position of the neighbor
  /// \p nghbrPos (w.r.t the boundary cell) and sets the ghost cell coordinates.
  ///
  /// \returns the local ghostCellId
  ///
  /// Precondition: cell at \p localBCellId has no neighbor in \nghbrPos
  CellIdx create_ghost_cell(const CellIdx bndryCellIdx, const SInd nghbrPos) {
    ASSERT(!is_valid(cells().neighbors(bndryCellIdx,nghbrPos)), "There is already a cell there!");

    /// 0) create a ghost cell without a globalBCellId
    const auto ghostCellIdx = cells().push_cell();
    cells().node_idx(ghostCellIdx) = invalid<NodeIdx>();

    /// 1) Set neighbor relationships:
    const auto oppositeNghbrPos = grid().nodes().opposite_nghbr_position(nghbrPos);
    cells().neighbors(bndryCellIdx,nghbrPos) = ghostCellIdx;
    cells().neighbors(ghostCellIdx,oppositeNghbrPos) = bndryCellIdx;

    /// 2) Set ghost cell position
    cells().x_center.row(ghostCellIdx) = ghost_cell_coordinates(ghostCellIdx);
    cells().length(ghostCellIdx) = cells().length(bndryCellIdx);
    return ghostCellIdx;
  }

  ///@}

  /// \name Initialization
  ///@{
  void init_solver_variables() {
    cfl_ = io::read<Num>(properties_,"CFL");
    dt_ = 0;
    time_ = 0;
    tEnd_ = io::read<Num>(properties_,"timeEnd");
    step_ = 0;
    forced_dt_ = 0;
    forced_dt_step_ = invalid<Ind>();
  }

  void create_local_cells() {
    // Create local cells:
    // create a local cell for each grid leaf cell
    // set the nodeIdx in the local cells
    // set the cIdx in the grid cells
    // set the cell coordinates
    {
      auto initialDomain = io::read<InitialDomain>(properties_,"initialDomain");
      for(auto nIdx : grid().nodes().leaf_nodes()) {
        auto xc = grid().cell_coordinates(nIdx);
        if(!initialDomain(xc)){ continue; }
        auto cIdx = cells().push_cell();
        cells().node_idx(cIdx) = nIdx;
        grid().nodes().cell_id(nIdx,solver_idx()) = cIdx;
        cells().x_center.row(cIdx) = xc;
        cells().length(cIdx) = grid().cell_length(nIdx);
      }
    }

    ASSERT(check_all_cells(),"solver cells / grid node links are wrong!");

    // set local nghbr ids: (localCellId,nIdx) -> globalNghbrIds -> localNghbrIds
    for(auto cIdx : cell_ids()) {
      auto nodeNghbrIds =  grid().nodes().template all_samelvl_nghbrs<strict>(node_idx(cIdx));
      for(auto nghbrPos : grid().nodes().nghbr_pos()) {
        auto nghbrIdx = nodeNghbrIds(nghbrPos);
        cells().neighbors(cIdx,nghbrPos)
            = is_valid(nghbrIdx) ? grid().nodes().cell_id(nghbrIdx,solver_idx())
                                 : invalid<CellIdx>();
      }
    }

    firstGC_ = CellIdx{cells().size()};

    create_ghost_cells();

    /// sort ghost cells by boundary id
    sort_gc();

    /// set distances to nghbrs
    {
      for(auto cIdx : cell_ids()) {
        for(const auto nghbrPos : grid().nodes().nghbr_pos()) {
          auto nghbrId = cells().neighbors(cIdx,nghbrPos);
          if(!is_valid(nghbrId)) { continue; } // \todo might not be necessary
          DBGV((cIdx)(nghbrId)(nghbrPos));
          cells().distances(cIdx,nghbrPos)
              = geometry::algorithm::distance(cells().x_center.row(cIdx).transpose(),
                                              cells().x_center.row(nghbrId).transpose());
          DBGV((cells().distances(cIdx,nghbrPos)));
        }
      }
    }
    ASSERT(check_all_nghbrs(),"internal cell nghbrIds do not agree with grid!");
  }

  /// Create Ghost Cells:
  ///
  /// \warning this only works for cutoff right now
  void create_ghost_cells() {
    auto noLeafCells = cells().size();
    for(auto boundaryCell : grid().boundary_cells(solver_idx())) {

      const auto bndryNodeIdx = boundaryCell.node_idx();

      const auto bndryCellIdx = grid().nodes().cell_id(bndryNodeIdx, solver_idx());
      ASSERT(is_valid(bndryCellIdx), "invalid bndryCellIdx!");
      ASSERT(node_idx(bndryCellIdx) == bndryNodeIdx, "solver and grid are not synchronized");
      ASSERT(is_valid(node_idx(bndryCellIdx)), "the global id has to be valid!");
      auto noCutBoundaries = boundaryCell.boundaries().size();

      switch(noCutBoundaries) {
        case 1: { // cell is only cut by one boundary
          const auto noMissingNghbrs = (cells().neighbors.row(bndryCellIdx).array()
                                        == invalid<CellIdx>()).count();
          // find in which direction the cell has no neighbor
          if(noMissingNghbrs == 0) {
            continue; // \todo remove at cut-off cleanup
          } else if(noMissingNghbrs == 1) {
            ASSERT( (cells().neighbors.row(bndryCellIdx).array()
                     == invalid<CellIdx>()).count() == 1,
                    "cell: localCellId = " << bndryCellIdx
                    << " (nodeIdx = " << node_idx(bndryCellIdx)
                    << ") is not missing one neighbor!");
            const CellIdxA<2*nd> nghbrs = cells().neighbors.row(bndryCellIdx);
            const auto nghbrIt = boost::find(nghbrs, invalid<CellIdx>());
            ASSERT(nghbrIt != boost::end(nghbrs), "missing neighbor not found!");
            const auto nghbrPos = nghbrIt - boost::begin(nghbrs);
            const auto ghostCellId = create_ghost_cell(bndryCellIdx,nghbrPos); // also sets coordinates!
            lastGC_ = ghostCellId;
            /// 2) Set boundary condition in the ghostcell
            cells().bc_idx(ghostCellId) = boundaryCell.boundaries()[0];
            break;
          } else { // if( noMissingNghbrs > 1 ) {
            /// \todo remove this cut-off hack
            const CellIdxA<2*nd> nghbrs = cells().neighbors.row(bndryCellIdx);
            auto nghbrIt = boost::find(nghbrs, invalid<CellIdx>());
            while(nghbrIt != boost::end(nghbrs)) {
              const auto nghbrPos = nghbrIt - boost::begin(nghbrs);
              const auto ghostCellId = create_ghost_cell(bndryCellIdx, nghbrPos); // also sets coordinates!
              lastGC_ = ghostCellId;
              cells().bc_idx(ghostCellId) = boundaryCell.boundaries()[0];
              nghbrIt = std::find(nghbrIt + 1, boost::end(nghbrs), invalid<CellIdx>());
            }
            break;
          }
        }
        case 2:
        case 3:
        case 4:
          /// \todo enable 5-6 only in 3D!
        case 5:
        case 6: { // cell is cut by more than one boundary ([2,6]) in 3D,
          // find where ghost cells should be added:
          std::vector<SInd> missingNghbrPos;
          for(const auto nghbrPos : grid().nodes().nghbr_pos()) {
            if(cells().neighbors(bndryCellIdx,nghbrPos) == invalid<CellIdx>()) {
              missingNghbrPos.push_back(nghbrPos);
            }
          }

          /// create ghost cells:
          std::vector<CellIdx> ghostCellIds;
          for(auto nghbrPos : missingNghbrPos) {
            ghostCellIds.push_back(create_ghost_cell(bndryCellIdx,nghbrPos));
          }

          /// each ghost cell has a negative signed distance for only one boundary
          std::vector<SInd> ghostCellBCIds;
          ghostCellBCIds.resize(ghostCellIds.size());
          for(auto ghostCellId : ghostCellIds) {
            auto ghostCellBoundaryId = invalid<SInd>();
            const typename Grid::CellVertices ghostCellVertices = cell_vertices(ghostCellId);

            for(auto gridBoundaryId : boundaryCell.boundaries()) {
              if(algorithm::all_of(ghostCellVertices(),[&](const NumA<nd>& x) {
                    return grid().boundaries()[gridBoundaryId].signed_distance(x) < 0;
                  })) {
                ASSERT(!is_valid(ghostCellBoundaryId),"overwritting boudnary id!");
                ghostCellBoundaryId = gridBoundaryId;
              }
            }
            ASSERT(is_valid(ghostCellBoundaryId),"boundary id not set!");
            cells().bc_idx(ghostCellId) = ghostCellBoundaryId;
          }
          break;
        }
        default: {
          TERMINATE(std::string("boundary cell is cut by ")
                    + std::to_string(noCutBoundaries)
                    + " boundaries!!");
        }
      }
    }
    auto newTotalNoCells = cells().size();
    std::cerr << "fv container | #of leafs: " << noLeafCells
              << " | #of ghosts: " << newTotalNoCells - noLeafCells
              << " | #of cells: " << newTotalNoCells << "\n";
  }

  /// \brief Sets the initial condition on the lhs
  void set_initial_condition() {
    auto initialCondition = io::read<InitialCondition>(properties_,"initialCondition");
    for(auto cIdx : internal_cells()) {
      const auto length = grid().cell_length(node_idx(cIdx));
      const auto x_c = cells().x_center.row(cIdx).transpose();
      const auto cellVol = geometry::cell::cartesian::volume<nd>(length);
      const auto cellInt = quadrature::integrate<nd>(initialCondition,x_c,length);
      const auto average = cellInt / cellVol;

      cells().lhs.row(cIdx) = average;
    }
  }
  ///@}
};

}} // namespace solver::fv

////////////////////////////////////////////////////////////////////////////////
} // hom3 namespace
////////////////////////////////////////////////////////////////////////////////
#undef ENABLE_DBG_
////////////////////////////////////////////////////////////////////////////////
#endif
