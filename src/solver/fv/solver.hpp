#ifndef HOM3_SOLVERS_FV_SOLVER_HPP_
#define HOM3_SOLVERS_FV_SOLVER_HPP_
////////////////////////////////////////////////////////////////////////////////
/// Includes:
#include "../../grid/grid.hpp"
#include "container.hpp"
#include "../../quadrature/quadrature.hpp"
/// Options:
#define ENABLE_DBG_ 0
#include "../../misc/dbg.hpp"
////////////////////////////////////////////////////////////////////////////////
/// Finite Volume cell:

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
   using NumFlux = NumFlux_;
  using Physics           = PhysicsTT<nd,Solver<nd,PhysicsTT,NumFlux/*,NumFluxTT*/>>;
  using physics_type      = typename Physics::type;
  static const SInd nvars = Physics::nvars;


  //  using NumFlux           = NumFluxTT<nd,Solver<nd,PhysicsTT/*,NumFluxTT*/>>;

  using CellContainer     = Container<nd,nvars>;
  using BoundaryCondition = std::function<void(Ind)>;
  using InitialCondition  = std::function<NumA<nvars>(const NumA<nd>)>;
  using InitialDomain     = std::function<bool(const NumA<nd>)>;

  using rhs               = rhs_tag;
  using lhs               = lhs_tag;
  ///@}

  // template<class T>
  // inline NumA<nvars> compute_num_flux(const Ind leftId, const Ind rightId,
  //                                     const SInd d, const Num dx, const Num dt) const {
  //   //return physics()->template ausm<T>(leftId,rightId,d,dx,dt);
  //   return physics()->template compute_num(leftId,rightId,d,dx,dt);
  //   //
  //   // return ausm<T>(leftId,rightId,d,dx,dt);
  // }

  /// \name Solver interface
  ///@{

  /// \brief Solver constructor
  ///
  /// Required properties are:
  /// - grid
  /// - maxNoCells
  /// - any extra properties required by the Physics class
  Solver(SInd solverId, io::Properties input) :
      Physics(input),
      solverId_(solverId),
      properties_(input),
      grid_(*io::read<Grid<nd>*>(input,"grid")),
      cells_(io::read<Ind>(input,"maxNoCells"))
      {}

  /// \brief Returns the solver id
  /// \complexity O(1)
  inline SInd solver_id() const { return solverId_; }
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
  const Ind& global_id(const Ind localId) const { return cells().globalId(localId); }
  /// \brief Returns the solver type-name (string)
  /// \warning for I\O only! For boolean operations please use solver_type() (unimplemented?)
  /// \complexity O(1)
  static std::string solver_type_name() { return "Fv"; }

  /// \todo splitt into multiple functions
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
  inline const Grid<nd>& grid() const { return grid_; }
  inline       Grid<nd>& grid()       { return grid_; }

  /// \brief Cell accessors
  inline       CellContainer& cells()       { return cells_; }
  inline const CellContainer& cells() const { return cells_; }

  /// \brief Conservative variables
  template<class T> inline NumA<nvars> Q(const Ind lId) const {
    return vars_(T(),lId);
  }
  template<class T> inline Eigen::Block<NumM<nvars>,1,nvars> Q(const Ind lId) {
    return vars_(T(),lId);
  }
  template<class T> inline Num  Q(const Ind lId, const SInd varId) const {
    return vars_(T(), lId, varId);
  }
  template<class T> inline Num& Q(const Ind lId, const SInd varId)       {
    return vars_(T(), lId, varId);
  }

private:
  /// \todo C&P code: REFACTOR: see grid/container.hpp (move to grid/range.hpp?)
  /// \brief Returns [RangeFilter] of existing node ids

  inline RangeFilter<Ind> valid() const {
    return RangeFilter<Ind>{[&](const Ind id){ return is_valid(id); }};
  }

  inline RangeFilter<Ind> ghosts() const {
    return RangeFilter<Ind>{[&](const Ind localId){ return is_ghost_cell(localId); }};
  }

  inline RangeFilter<Ind> not_ghosts() const {
    return RangeFilter<Ind>{[&](const Ind localId){ return !is_ghost_cell(localId); }};
  }

  inline RangeFilter<Ind> boundary() const {
    return RangeFilter<Ind>{[&](const Ind localId){ return is_valid(cells().bcId(localId)); }};
  }

  inline RangeFilter<Ind> boundary(const SInd bcId) const {
    return RangeFilter<Ind>{[&](const Ind localId){ return cells().bcId(localId) == bcId; }};
  }

  inline RangeTransformer<Ind,Ind> local_to_global() const {
    return {[&](const Ind localId){ return global_id(localId); }};
  }

  /// \brief Range of all local ids (previously all_cells)
  /// \warning This includes local ids without a global id
  inline auto local_ids() const -> Range<Ind> {
    return cells().all_cells();
  }

  inline auto ghost_cells() const -> FRange<Ind> {
      return local_ids() | ghosts();
  }

  inline auto boundary_cells() const -> FRange<Ind> {
      return local_ids() | boundary();
  }

  inline auto boundary_cells(const SInd bcId) const -> FRange<Ind> {
      return local_ids() | boundary(bcId);
  }

public:

  /// \brief Range of all internal cells
  /// i.e. all cells which are not ghost cells
  inline auto internal_cells() const -> FRange<Ind> {
      return local_ids() | not_ghosts();
  }

  /// \brief Range of global ids
  ///
  /// \todo This is not 100% right: local_ids should be maped using
  /// local_to_global into global_ids and then filtered by validity
  inline auto global_ids() const -> decltype(internal_cells() | local_to_global()) {
    return internal_cells() | local_to_global();
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
  /// Note: ghost cells can only have one boundary cell associated with them,
  /// right?
  inline auto boundary_cell_id(const Ind ghostCellId) const -> std::tuple<Ind,SInd> {
    for(SInd nghbrPos = 0; nghbrPos < no_nghbr_pos(); ++nghbrPos) {
      if(is_valid(cells().nghbrs(ghostCellId,nghbrPos))) {
        return std::make_tuple(cells().nghbrs(ghostCellId,nghbrPos),nghbrPos);
      }
    }
    /// this actually doesnt check for multiple nghbrs which can happen...
    ASSERT(false,"ghostcell are required to have only one neighbor!");
    return std::make_tuple(iInd(),iSInd());
  }

  ///@}


  /// \todo remove!
  void boundary_condition(SInd , Ind ) const { }

  /// \todo Have a general function, and enable_if a physics specific function
  /// per physics class (maybe write this function in the phyiscs class itself?)
  /// \name Output
  ///@{
  friend void write_domain(Solver* solver) {
    std::string fName = solver->solver_type_name() + "_"
                        + solver->physics_name() + "_"
                        + "Id" + std::to_string(solver->solver_id())
                        + "_" + std::to_string(solver->step());
    solver->apply_bcs();
    write_domain(fName,solver);
  }

  friend void write_domain(const std::string fName, const Solver* solver) {

    io::Vtk<nd,io::format::binary> out(
        io::StreamableDomain<nd>( // the following is not necessary but really cool:
            [&](){return AnyRange<Ind>(algorithm::join(solver->local_ids() | solver->not_ghosts(),
                                                       solver->local_ids() | solver->ghosts())); },
            [&](const Ind lId){ return solver->cell_vertices(lId);}),
        fName, io::precision::standard());

    out << io::stream("locallCellIds",1,[](const Ind cId, const SInd){ return cId; });
    out << io::stream("local_nghbrs",solver->grid().nodes().no_samelvl_nghbr_pos(),
                      [&](const Ind lId, const SInd pos){
                        return solver->cells().nghbrs(lId,pos);
                      });
    out << io::stream("ghost_cell",1,[&](const Ind lId, const SInd){
        return solver->is_ghost_cell(lId) ? 1 : iInd();
      });
    out << io::stream("bcId",1,[&](const Ind lId, const SInd) -> Ind {
        return solver->cells().bcId(lId);
      });

    out << io::stream(Solver::V::cv_names,nvars,[&](const Ind cId, const SInd d){
        return solver->cells().vars(cId, d);
      });

    solver->physics()->template physics_output(out);
  }
  ///@}

private:

  /// \name General data
  ///@{

  /// Stores the solver id
  const SInd solverId_;
  /// Solver properties
  io::Properties properties_;
  /// Solver grid
  Grid<nd>& grid_;
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
  Ind firstGC_;

  /// Offset to the last ghost cell in the container
  Ind lastGC_;

  // Avoids insanity while working with Physics CRTP
  friend PhysicsTT<nd,Solver<nd,PhysicsTT,NumFlux/*,NumFluxTT*/>>;

  /// \brief Access to Physics class
  // Note: think of it as the inverse of a CRTP
  inline       Physics* physics()       { return static_cast<Physics*>(this); }
  inline const Physics* physics() const { return static_cast<const Physics*>(this); }

  /// \name Tag-dispatching for conservative variables (rhs/lhs)
  ///@{
  inline const Num& vars_(rhs, const Ind lId, const SInd varId) const {
    return cells().rhs(lId,varId);
  }
  inline       Num& vars_(rhs, const Ind lId, const SInd varId)       {
    return cells().rhs(lId,varId);
  }
  inline NumA<nvars> vars_(rhs,const Ind lId) const {
    return cells().rhs().row(lId);
  }
  inline Eigen::Block<NumM<nvars>,1,nvars> vars_(rhs,const Ind lId) {
    return cells().rhs().row(lId);
  }
  inline const Num& vars_(lhs, const Ind lId, const SInd varId) const {
    return cells().vars(lId,varId);
  }
  inline       Num& vars_(lhs, const Ind lId, const SInd varId)       {
    return cells().vars(lId,varId);
  }
  inline NumA<nvars> vars_(lhs, const Ind lId) const {
    return cells().vars().row(lId);
  }
  inline Eigen::Block<NumM<nvars>,1,nvars> vars_(lhs, const Ind lId) {
    return cells().vars().row(lId);
  }
  ///@}

  /// \name Numerical functions
  ///@{

  /// \brief Computes the numerical flux
  /// \todo Numerical flux function should be configurable (using a TP?)
  template<class T>
  inline NumA<nvars> num_flux(const Ind localId) const {
    DBG("first term in rhs:");
    DBGV((localId));
    NumA<nvars> result = NumA<nvars>::Zero();
    const Num dx = cells().length(localId);
    for(SInd d = 0; d < nd; ++d) {
      const SInd nghbrM = d*2;
      const SInd nghbrP = nghbrM + 1;
      const Ind nghbrMId = cells().nghbrs(localId,nghbrM);
      const Ind nghbrPId = cells().nghbrs(localId,nghbrP);
      const auto flux_m = physics()->template compute_num_flux<T,NumFlux>(nghbrMId,localId,d,dx,dt());
      const auto flux_p = physics()->template compute_num_flux<T,NumFlux>(localId,nghbrPId,d,dx,dt());
      result += dt() / dx * (flux_m - flux_p);
      DBGV((localId)(d)(dx)(dt())(nghbrMId)(nghbrPId)(flux_m)(flux_p)(result)
           (Q<T>(nghbrMId))(Q<T>(localId))(Q<T>(nghbrPId)));
    }
    return result;
  }

  /// \brief Performs an Euler-Forward step for all cells in range \p cells
  template<class CellRange>
  inline void ef(CellRange&& cells) {
    apply_bcs();
    for(const Ind cId : cells) {
      Q<rhs>(cId) = Q<lhs>(cId) + num_flux<lhs>(cId).transpose();
    }
    for(auto cId : cells) {
      Q<lhs>(cId) = Q<rhs>(cId);
    }
  }

  /// \brief Performs a RK2 step for all cells in range \p cells
  ///
  /// \todo apply_bcs should be parametrized to with rhs/lhs to decide when what
  /// should be applied
  template<class CellRange>
  inline void rk2(CellRange&& cells) {
    apply_bcs();
    for(auto cId : cells) {
      Q<rhs>(cId) = Q<lhs>(cId) + num_flux<lhs>(cId).transpose();
    }

    apply_bcs();
    for(auto cId : cells) {
      Q<lhs>(cId) = 0.5 * (Q<lhs>(cId) + Q<rhs>(cId)
                           + num_flux<rhs>(cId).transpose());
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
    Num dt = std::numeric_limits<Num>::max();
    for(const Ind cId : internal_cells())  {
      dt = std::min(dt,physics()->template compute_dt<lhs>(cId));
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
    SInd bcId
        = cells().bcId(container::sequential::algorithm::find_if(0,cells().size(),[&](Ind i){
              return cells().bcId(i) != iSInd(); }));
    auto eql_bcId = [&](Ind i) { return cells().bcId(i) == bcId; };
    Ind firstGhostCell = container::sequential::algorithm::find_if(cells(),eql_bcId);
    Ind lastGhostCell = firstGhostCell;
    //DBGV_ON((bcId)(firstGhostCell));
    for(auto& boundary: grid().boundaries(solver_id())) {
      ++bcId;
      firstGhostCell = lastGhostCell;
      lastGhostCell = container::sequential::algorithm::find_if(firstGhostCell,cells().size(),eql_bcId);
      //DBGV_ON((bcId)(firstGhostCell)(lastGhostCell));
      auto GCRange = Range<Ind>{firstGhostCell,lastGhostCell};
      // std::cerr << physics()->physics_name() << " | bcId: " << bcId
      //           << " | range: [" <<  firstGhostCell << ", " << lastGhostCell << ")"
      //           << std::endl;
      boundary.condition(solver_id())(GCRange);
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

  /// \brief Checks that all global nghbrs of local cell \p localId agree with
  /// those in the grid
  void check_nghbrs(const Ind localId) const {
    const auto globalId = global_id(localId);
    ASSERT(is_valid(globalId), "you can only check cells that have valid global ids!");
    const auto globalNghbrIds = grid().nodes().template all_samelvl_nghbrs<strict>(globalId);
    SInd nghbrPos = 0;
    for(auto globalNghbrId : globalNghbrIds) {
      if(is_valid(globalNghbrId) // is there a global neighbor there belonging to the same domain?
         && grid().nodes().has_container(globalNghbrId,solver_id())) {
        const auto localNghbrId = grid().local_id(globalNghbrId,solver_id());
        ASSERT(cells().nghbrs(localId,nghbrPos) == localNghbrId,
               "cell: localId = " << localId << " (globalId = " << global_id(localId) <<  ") has a wrong nghbr in pos " << nghbrPos
               << "(gridId: " << localNghbrId << ", local solverId: "
               << cells().nghbrs(localId,nghbrPos));
      }
      ++nghbrPos;
    }
  }

  /// \brief Checks neighbros for all internal cells
  bool check_all_nghbrs() const {
    for(auto localId : internal_cells()) {
      check_nghbrs(localId);
    }
    return true;
  }

  //@}

  /// \name Grid/Spatial functions (some of it is for convenience)
  /// \todo Refactor this functionality into non-member functions of grid?
  ///@{

  /// \brief Number of samelvl neighbor positions
  constexpr SInd no_nghbr_pos() const { return grid().nodes().no_samelvl_nghbr_pos(); }

  /// \brief Number of local cell neighbors (includes ghost-cells)
  /// \todo recycle somehow the hiearchical container algorithm?
  SInd no_nghbrs(const Ind cellId) const {
    SInd noNghbrs = 0;
    for(SInd nghbrPos = 0; nghbrPos < no_nghbr_pos(); ++nghbrPos) {
      if(is_valid(cells().nghbrs(cellId,nghbrPos))) {
        ++noNghbrs;
      }
    }
    return noNghbrs;
  }

  /// \brief Gets the cell length from the grid
  /// \todo refactor into: set_length(cId) method
  /// \todo cell_length should return the length stored in the container O(1)
  /// instead of O(logN)
  inline Num cell_length(const Ind localCellId) const {
    ASSERT(is_valid(localCellId),"works only for cells with global id!");
    return grid().cell_length(global_id(localCellId));
  }

  /// \brief Lazy range of all cell vertices
  /// \todo Check, I think the range is strict now
  auto cell_vertices(const Ind localCellId) const -> grid::CellVertices<nd> {
    const auto x_cell = cells().xc().row(localCellId);

    const Num length = !is_ghost_cell(localCellId)
      ? cell_length(localCellId) // cell is a global cell -> get its length
      : [&](){ // cell is a ghost cell -> get its boundary cell's length
      Ind bndryCellId = iInd();
      std::tie(bndryCellId,std::ignore) = boundary_cell_id(localCellId);
      return cell_length(bndryCellId);
    }();

    return { localCellId, grid().cell_vertices_coords(length, x_cell) };
  }
  ///@}

  /// \name Ghost-cell related
  ///@{

  inline bool is_ghost_cell(const Ind localId) const {
    return !is_valid(global_id(localId));
  }

  /// \brief Computes the cell-center coordinates of the ghost cell \p
  /// ghostCellId
  NumA<nd> ghostcell_coords(const Ind ghostCellId)  {
    ASSERT(no_nghbrs(ghostCellId) == 1, "Invalid ghost cell!");

    /// Find local boundary cell id and the ghost cell position w.r.t. the
    /// boundary cell
    Ind localBndryId = iInd();
    SInd position_wrt_ghostCell = iSInd();
    std::tie(localBndryId,position_wrt_ghostCell) = boundary_cell_id(ghostCellId);
    ASSERT(is_valid(position_wrt_ghostCell), "ERROR!");
    ASSERT(is_valid(localBndryId),"ERROR!");

    /// \todo: opposite nghbr position should be a non-member function!
    const SInd position_wrt_boundary = grid().nodes().opposite_nghbr_position(position_wrt_ghostCell);
    ASSERT(is_valid(position_wrt_boundary),"ERROR!");

    auto length = cell_length(localBndryId);
    NumA<nd> x_gc = cells().xc().row(localBndryId);
    x_gc += length * grid().nghbr_rel_pos(position_wrt_boundary).template cast<Num>();
    return x_gc;
  }

  /// \brief Sorts ghost cells by boundary condition id and returns the id of
  /// the first ghost cell
  Ind sort_gc() {

    /// Sort ghost cells by bcId:
    auto valid_bcId = [&](Ind i) { return is_valid(cells().bcId(i)); };
    const Ind firstGhostCell = container::sequential::algorithm::find_if(cells(),valid_bcId);

    auto sort_bcIds = [&](typename CellContainer::value_type a,
                          typename CellContainer::value_type b) {
      return a.bcId() < b.bcId();
    };
    std::sort(cells().begin()+firstGhostCell,cells().end(),sort_bcIds);

    /// Correct nghbrs in boundary cells:
    for(auto ghostId : ghost_cells()) {
      for(auto nghbrPos : grid().nodes().nghbr_pos()) {
        auto nghbrId = cells().nghbrs(ghostId,nghbrPos);
        if(is_valid(nghbrId)) {
          auto oppositeNghbrPos = grid().nodes().opposite_nghbr_position(nghbrPos);
          cells().nghbrs(nghbrId,oppositeNghbrPos) = ghostId;
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
  Ind create_ghost_cell(const Ind localBCellId, const SInd nghbrPos) {
    ASSERT(!is_valid(cells().nghbrs(localBCellId,nghbrPos)), "There is already a cell there!");

    /// 0) create a ghost cell without a globalBCellId
    cells().push_cell();
    const Ind ghostCellId = cells().back();
    global_id_(ghostCellId) = iInd();

    /// 1) Set neighbor relationships:
    const SInd oppositeNghbrPos = grid().nodes().opposite_nghbr_position(nghbrPos);
    cells().nghbrs(localBCellId,nghbrPos) = ghostCellId;
    cells().nghbrs(ghostCellId,oppositeNghbrPos) = localBCellId;

    /// 2) Set ghost cell position
    cells().xc().row(ghostCellId) = ghostcell_coords(ghostCellId);
    cells().length(ghostCellId) = cells().length(localBCellId);
    return ghostCellId;
  }

  ///@}


  /// \name Cell container manipulation
  ///@{

  /// Maps local ids with global ids (returns a writable reference!)
  Ind& global_id_(const Ind localId)       { return cells().globalId(localId); }

  /// \brief Sets cell-center coordinates of cell \p localCellId
  void set_cell_coordinates(const Ind localCellId, const NumA<nd> xc) {
    ASSERT(is_valid(global_id(localCellId)), "doesn't work with non-global cells!");
    for(auto d : grid().dims()) {
      cells().xc(localCellId,d) = xc(d);
    }
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
    forced_dt_step_ = iInd();
  }

  void create_local_cells() {
    // Create local cells:
    // create a local cell for each grid leaf cell
    // set the globalIds in the local cells
    // set the localIds in the grid cells
    // set the cell coordinates
    {
      auto initialDomain = io::read<InitialDomain>(properties_,"initialDomain");
      Ind localId = 0;
      for(auto globalCellId : grid().nodes().leaf_nodes()) {
        auto xc = grid().cell_coordinates(globalCellId);
        if(!initialDomain(xc)){ continue; }
        cells().push_cell();
        global_id_(localId) = globalCellId;
        grid().local_id(globalCellId,solver_id()) = localId;
        set_cell_coordinates(localId,xc);
        cells().length(localId) = grid().cell_length(globalCellId);
        ++localId;
      }
    }


    // set local nghbr ids: (localCellId,globalCellId) -> globalNghbrIds -> localNghbrIds
    for(auto localId : local_ids()) {
      auto globalId = global_id(localId);
      cells().nghbrs().row(localId)
          = grid().nodes().template all_samelvl_nghbrs<strict>(globalId).unaryExpr([&](const Ind i) {
              return is_valid(i) ? grid().local_id(i,solver_id()) : i;
            });
    }

    firstGC_ = cells().size();

    create_ghost_cells();

    /// sort ghost cells by boundary id
    sort_gc();

    // for(auto id : local_ids()) {
    //   auto bcId = cells().bcId(id);
    //   if(is_valid(bcId)) {
    //     std::cerr << "cId: " << id << " bcId: " << bcId << std::endl;
    //   }
    // }

    /// set distances to nghbrs
    {
      for(auto localId : local_ids()) {
        for(SInd nghbrPos = 0; nghbrPos < 2*nd; ++nghbrPos) {
          auto nghbrId = cells().nghbrs(localId,nghbrPos);
            if(!is_valid(nghbrId)) { continue; } // \todo might not be necessary
            DBGV((localId)(nghbrId)(nghbrPos));
            cells().dx(localId,nghbrPos)
                = geometry::algorithm::distance(cells().xc().row(localId).transpose(),
                                                cells().xc().row(nghbrId).transpose());
            DBGV((cells().dx(localId,nghbrPos)));
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
    for(auto boundaryCell : grid().boundary_cells(*this)) {

      const auto globalBCellId = boundaryCell.global_id();

      const auto localBCellId = grid().nodes().local_id(globalBCellId, solver_id());
      ASSERT(global_id(localBCellId) == globalBCellId, "solver and grid are not synchronized");
      ASSERT(is_valid(global_id(localBCellId)), "the global id has to be valid!");
      auto noCutBoundaries = boundaryCell.boundaries().size();

      switch(noCutBoundaries) {
        case 1: { // cell is only cut by one boundary
          auto noMissingNghbrs = (cells().nghbrs().row(localBCellId).array() == iInd()).count();
          // find in which direction the cell has no neighbor
          if(noMissingNghbrs == 0) {
            continue; // \todo remove at cut-off cleanup
          } else if(noMissingNghbrs == 1) {
            ASSERT( (cells().nghbrs().row(localBCellId).array() == iInd()).count() == 1,
                    "cell: localCellId = " << localBCellId
                    << " (globalId = " << global_id(localBCellId)
                    << ") is not missing one neighbor!");
            const IndA<2*nd> nghbrs  = cells().nghbrs().row(localBCellId);
            auto pos = boost::find(nghbrs,iInd());
            ASSERT(pos != boost::end(nghbrs), "missing neighbor not found!");
            auto nghbrPos = pos - boost::begin(nghbrs);
            auto ghostCellId = create_ghost_cell(localBCellId,nghbrPos); // also sets coordinates!
            lastGC_ = ghostCellId;
            /// 2) Set boundary condition in the ghostcell
            cells().bcId(ghostCellId) = boundaryCell.boundaries()[0];
            break;
          } else { // if( noMissingNghbrs > 1 ) {
            /// \todo remove this cut-off hack
            const IndA<2*nd> nghbrs  = cells().nghbrs().row(localBCellId);
            auto pos = boost::find(nghbrs,iInd());
            while(pos != boost::end(nghbrs)) {
              auto nghbrPos = pos - boost::begin(nghbrs);
              auto ghostCellId = create_ghost_cell(localBCellId,nghbrPos); // also sets coordinates!
              lastGC_ = ghostCellId;
              cells().bcId(ghostCellId) = boundaryCell.boundaries()[0];
              //pos = boost::find(nghbrs,iInd());
              pos = std::find(pos+1,boost::end(nghbrs),iInd());
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
          for(SInd nghbrPos = 0; nghbrPos < 2*nd; ++nghbrPos) {
            if(cells().nghbrs(localBCellId,nghbrPos) == iInd()) {
              missingNghbrPos.push_back(nghbrPos);
            }
          }

          /// create ghost cells:
          std::vector<Ind> ghostCellIds;
          for(auto nghbrPos : missingNghbrPos ) {
            ghostCellIds.push_back(create_ghost_cell(localBCellId,nghbrPos));
          }

          /// each ghost cell has a negative signed distance for only one boundary
          std::vector<SInd> ghostCellBCIds;
          ghostCellBCIds.resize(ghostCellIds.size());
          for(auto ghostCellId : ghostCellIds) {
            auto ghostCellBoundaryId = iSInd();
            const grid::CellVertices<nd> ghostCellVertices = cell_vertices(ghostCellId);

            for(auto gridBoundaryId : boundaryCell.boundaries()) {
              if(algorithm::all_of(ghostCellVertices(),[&](const NumA<nd>& x) {
                    return grid().boundaries()[gridBoundaryId].signed_distance(x) < 0;
                  })) {
                ASSERT(!is_valid(ghostCellBoundaryId),"overwritting boudnary id!");
                ghostCellBoundaryId = gridBoundaryId;
              }
            }
            ASSERT(is_valid(ghostCellBoundaryId),"boundary id not set!");
            cells().bcId(ghostCellId) = ghostCellBoundaryId;
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

  void set_initial_condition() {
    auto initialCondition = io::read<InitialCondition>(properties_,"initialCondition");
    for(auto localId : internal_cells()) {
      auto length = grid().cell_length(global_id(localId));
      const NumA<nd> x_cell = cells().xc().row(localId);
      auto cellVol = geometry::cell::cartesian::volume<nd>(length);
      auto cellInt = quadrature::integrate<nd>(initialCondition,x_cell,length);
      auto average = cellInt / cellVol;
      //auto cellVars = initialCondition(cells().xc().row(localId));
      //DBGV_ON((localId)(length)(cellVol)(cellInt)(average)(x_cell));
      for(SInd v = 0; v < nvars; ++v) {
        cells().vars(localId,v) = average(v); // gauss-q
        //cells().vars(localId,v) = cellVars(v); // simple
      }
    }
    //apply_bcs();
  }
  ///@}
};

}} // namespace solver::fv

////////////////////////////////////////////////////////////////////////////////
#undef ENABLE_DBG_
////////////////////////////////////////////////////////////////////////////////
#endif
