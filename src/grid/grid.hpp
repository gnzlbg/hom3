#ifndef HOM3_GRID_HPP_
#define HOM3_GRID_HPP_
////////////////////////////////////////////////////////////////////////////////
/// \file \brief Contains the grid data-structure
///
/// \todo Replace warnings by compile-time/run-time assertions
////////////////////////////////////////////////////////////////////////////////
/// Includes:
#include "../containers/hierarchical.hpp"
#include "generation.hpp"
#include "cellvertex.hpp"
#include "../geometry/boundary.hpp"
/// Options:
#define ENABLE_DBG_ 0
#include "../misc/dbg.hpp"
////////////////////////////////////////////////////////////////////////////////

namespace grid {

/// \name Grid tags
///@{

/// \name Constructor tags
///@{
/// \brief Initializes the grid on construction
struct initialize_t {};
static constexpr initialize_t initialize{};
///@}

/// \brief #of space dimensions tag
template<SInd nd> struct dim {};

///@}

/// \name child_rel_pos stencils
///@{
static constexpr std::array<SInt,4*2> child_rel_pos_2d_arr {{
    // x   y
      -1, -1, // pos: 0
       1, -1, // pos: 1
      -1,  1, // pos: 2
       1,  1  // pos: 3
    }};

static constexpr std::array<SInt,8*3> child_rel_pos_3d_arr {{
    // x   y   z
      -1, -1, -1 , // pos: 0
       1, -1, -1 , // pos: 1
      -1,  1, -1 , // pos: 2
       1,  1, -1 , // pos: 3
      -1, -1,  1 , // pos: 4
       1, -1,  1 , // pos: 5
      -1,  1,  1 , // pos: 6
       1,  1,  1   // pos: 7
    }};
///@}

/// \name nghbr_rel_pos stencils
///@{
static constexpr std::array<SInt,2*4> nghbr_rel_pos_2d_arr {{
    // x   y   z
      -1, 0, // pos: 0
       1, 0, // pos: 1
       0,-1, // pos: 2
       0, 1  // pos: 3
    }};

static constexpr std::array<SInt,3*6> nghbr_rel_pos_3d_arr {{
    // x   y   z
       -1, 0, 0, // pos: 0
        1, 0, 0, // pos: 1
        0,-1, 0, // pos: 2
        0, 1, 0, // pos: 3
        0, 0,-1, // pos: 4
        0, 0, 1  // pos: 5
    }};
///@}

/// \name vertex_rel_pos stencils
///@{
static constexpr std::array<SInt,4*2> vertex_pos_2d_arr {{
    // x   y   z
      -1, -1 , // pos: 0
       1, -1 , // pos: 1
       1,  1 , // pos: 2
      -1,  1   // pos: 3
    }};

static constexpr std::array<SInt,8*3> vertex_pos_3d_arr {{
    // x   y   z
      -1, -1, -1 , // pos: 0
       1, -1, -1 , // pos: 1
       1,  1, -1 , // pos: 2
      -1,  1, -1 , // pos: 3
      -1, -1,  1 , // pos: 4
       1, -1,  1 , // pos: 5
       1,  1,  1 , // pos: 6
      -1,  1,  1   // pos: 7
    }};
///@}

template<int nd> struct RootCell {
  RootCell(const RootCell<nd>&) = default;
  RootCell(const NumA<nd>& min, const NumA<nd>& max) : length(math::eps) {
    for(Ind d = 0; d != nd; ++d) {
      const Num lengthD = max(d) - min(d);
      ASSERT(lengthD > math::eps, "Negative length not allowed!");
      length = std::max(length, lengthD);
      coordinates(d) = min(d) + 0.5 * lengthD;
    }
  }

  static const int nDim = nd;
  Num length;
  NumA<nd> coordinates;
  NumA<nd> x_min() const {
    NumA<nd> m;
    for(SInd d = 0; d < nd; ++d) {
      m(d) = coordinates(d) - 0.5 * length;
    }
    return m;
  }
  NumA<nd> x_max() const {
    NumA<nd> m;
    for(SInd d = 0; d < nd; ++d) {
      m(d) = coordinates(d) + 0.5 * length;
    }
    return m;
  }
};

} // grid namespace

/// \brief Hierarchical Cartesian Grid data-structure
/// - couples the spatial information (e.g. coordinates) with the connectivity
/// graph (e.g. neighbor/parent-child relationships)
/// - manages the boundaries/geometry of each sub domain
///
/// \todo rename to CartesianHSP
template <SInd nd> struct CartesianHSP {

  using This = CartesianHSP<nd>;
  using Boundary = boundary::Interface<nd>;
  using Boundaries = std::vector<Boundary>;
  using dim = grid::dim<nd>;
  using D2 = grid::dim<2>;
  using D3 = grid::dim<3>;
  using Connectivity = typename container::Hierarchical<nd>;
  /// \name Construction
  ///@{

  /// \brief Construct a grid from a set of input properties
  CartesianHSP(io::Properties input)
      : properties_(input),
        nodes_(io::read<Ind>(properties_,"maxNoGridNodes"),
               io::read<SInd>(properties_,"maxNoContainers")),
        rootCell_(io::read<grid::RootCell<nd>>(properties_,"rootCell")),
        ready_(false)
  { TRACE_IN_(); TRACE_OUT(); }

  /// \brief Constructs and initializes a grid from a set of input properties
  CartesianHSP(io::Properties input, grid::initialize_t) : CartesianHSP(input) {
    TRACE_IN_();
    initialize();
    TRACE_OUT();
  }

  /// \brief Initializes the grid, i.e. reads the boundaries and the mesh
  /// generator from the properties if they haven't been set yet and generates
  /// the mesh
  void initialize() {
    if(boundaries().empty()) {
      read_boundaries();
    }
    if(!meshGeneration_) {
      read_mesh_generator();
      generate_mesh();
    }
    ready_ = true;
  }

  /// \brief Reads boundaries from input set
  ///
  /// \warning Overwrites previously-set boundaries!
  void read_boundaries() { read_boundaries(properties_); }
  void read_boundaries(const io::Properties& boundaryProperties) {
    io::read(boundaryProperties,"boundaries",boundaries_);
  }

  /// \brief Reads mesh-generator from input set
  ///
  /// \warning Overwrites previously-set mesh-generator!
  void read_mesh_generator() {
    io::read(properties_,"meshGeneration",meshGeneration_);
  }

  /// \brief Sets functor \p meshGenerator as the grid mesh generator
  ///
  /// \warning Overwrites previously-set mesh-generator!
  template<class MeshGenerator>
  void set_mesh_generator(MeshGenerator&& meshGenerator) {
    meshGeneration_ = meshGenerator;
  }

  /// \brief Generates a mesh using the grid's mesh-generator
  ///
  /// \warning Doesn't check if a mesh-generator has been set
  /// \warning Doesn't check if a grid has been generated already!
  void generate_mesh() {
    TRACE_IN_();
       meshGeneration_(*this);
    TRACE_OUT();
  };

  ///@}

  /// \name Access
  ///@{

  /// \brief Range over space dimensions
  inline auto dims() const -> Range<SInd> { return {SInd(0),nd}; }

  /// \brief Mesh topology
  inline const Connectivity& nodes() const { return nodes_; }
  inline       Connectivity& nodes()       { return nodes_; }

  ///@}

  /// \name Spatial information
  ///@{

  /// \brief #of spatial dimensions
  static constexpr SInd no_dimensions() {
    static_assert(nd == 2 || nd == 3,
                  "Error: unsupported #of spatial dimensions!");
    return nd;
  }
  /// \brief #of verticles per cell
  static constexpr SInd no_edge_vertices() { return math::ct::ipow(2u,nd); }

  /// \brief Relative position (-1 or +1) of child located at \p childPos w.r.t
  /// its parent along the \p d -axis.
  ///
  /// That is:
  ///              __________________________
  ///            /|   pos: 6   |   pos: 7  /|
  ///           / | (-1,+1,+1) | (+1,+1,+1) |
  ///          /  |____________|____________|
  ///         /   |   pos: 4   |   pos: 5   |
  ///        /    | (-1,-1,+1) | (+1,-1,+1) |
  ///       /     |____________|____________|
  ///      /     /                   /     /
  ///     /_____/___________________/     /
  ///    |   pos: 2   |   pos: 3   |     /    y (1) ^
  ///    | (-1,+1,-1) | (+1,+1,-1) |    /           |     ^ z (2)
  ///    |____________|____________|   /            |    /
  ///    |   pos: 0   |   pos: 1   |  /             |  /
  ///    | (-1,-1,-1) | (+1,-1,-1) | /              |/
  ///    |____________|____________|/               o-------> x (0)
  ///
  ///
  /// \todo some sort of constexpr assertion for valid childPos and d
  ///
  static constexpr SInt child_rel_pos(const SInd childPos, const SInd d) {
    return child_rel_pos_(childPos,d,dim());
  }

  /// \brief Relative position of child localted at \p childPos w.r.t. its parent
  ///
  /// \returns vector of relative positions
  static constexpr SIntA<nd> child_rel_pos(const SInd childPos) {
    return child_rel_pos_(childPos,dim());
  }

  /// \brief Relative position of neighbor in position \p nghbrPos w.r.t. cell
  /// along the \p -axis.
  ///
  /// That is:
  ///
  ///                      _____________
  ///                     |   pos: 5   |
  ///                     | ( 0, 0,+1) |
  ///                     |____________|
  ///                  __/__________   /
  ///                 |   pos: 3   |  /
  ///                 | ( 0,+1, 0) | /
  ///     ____________|____________|/____________
  ///    |   pos: 0   |            |   pos: 1   |  y (1) ^
  ///    | (-1, 0, 0) |    cell    | ( 1, 0, 0) |        |     ^ z (2)
  ///    |____________|____________|____________|        |    /
  ///                /|   pos: 2  /|                     |  /
  ///              /  | ( 0,-1, 0) |                     |/
  ///            /    |____________|                     o-------> x (0)
  ///          /            /
  ///        /_____________/
  ///        |   pos: 4   |
  ///        | ( 0, 0,-1) |
  ///        |____________|
  ///
  ///
  static constexpr SInt nghbr_rel_pos(const SInd nghbrPos, const SInd d) {
    return ghbr_rel_pos_(nghbrPos,d,dim());
  }

  /// \brief Relative position of neighbor localted at \p nghbrPos w.r.t. its parent
  ///
  /// \returns vector of relative positions
  static constexpr SIntA<nd> nghbr_rel_pos(const SInd nghbrPos) {
    return nghbr_rel_pos_(nghbrPos,dim());
  }

  /// \brief Relative position of vertex in position \p vertexId (sorted in
  /// counter clock-wise order) w.r.t. to cell along the \p -axis.
  ///
  /// That is:
  ///              o________________________o
  ///            /|   pos: 7   |   pos: 6  /|
  ///           / | (-1,+1,+1) | (+1,+1,+1) |
  ///          /  |____________|____________|
  ///         /   |   pos: 4   |   pos: 5   |
  ///        /    | (-1,-1,+1) | (+1,-1,+1) |
  ///       /     o____________|____________o
  ///      /     /                   /     /
  ///     o_____/__________________o/     /
  ///    |   pos: 3   |   pos: 2   |     /    y (1) ^
  ///    | (-1,+1,-1) | (+1,+1,-1) |    /           |     ^ z (2)
  ///    |____________|____________|   /            |    /
  ///    |   pos: 0   |   pos: 1   |  /             |  /
  ///    | (-1,-1,-1) | (+1,-1,-1) | /              |/
  ///    o____________|____________o/               o-------> x (0)
  ///
  /// Note:
  ///  - the positions correspond to the vertices.
  ///  - the stencil is the same one as that for the relative child position,
  ///    just sorted in a different order.
  static constexpr SInt vertex_pos(const NodeIdx childIdx, const SInt d) {
    return vertex_pos_(childIdx,d,dim());
  }

  /// \brief Length of cells at \p level
  Num cell_length_at_level(const Ind level) const {
    return rootCell_.length / std::pow(2,level);
  }

  /// \brief Length of cell at node \p nIdx
  Num cell_length(const NodeIdx nIdx) const {
    return cell_length_at_level(nodes().level(nIdx));
  }

  /// \brief Centroid coordinates of cell at node \p nIdx
  ///
  /// Note: the coordinates are computed from the root cell's coordinates.
  /// \complexity O(L) - linear in the number of levels.
  /// \complexity O(logN) - logarithmic in the number of cells.
  NumA<nd> cell_coordinates(const NodeIdx nIdx) const {
    TRACE_IN((nIdx));
    const Ind level = nodes().level(nIdx);
    if( level != 0 ) {
      const auto pIdx      = nodes().parent(nIdx);
      const auto posInParent  = nodes().position_in_parent(nIdx);
      const auto relLength = 0.25 * cell_length(pIdx);
      const NumA<nd> relativePosition
          = cell_coordinates(pIdx)
          + child_rel_pos(posInParent).template cast<Num>() * relLength;
      TRACE_OUT();
      return relativePosition;
    } else {
      TRACE_OUT();
      return rootCell_.coordinates;
    }
  }

  /// \todo rename cell_id to node_id
  struct CellVertices {
    static constexpr SInd no_vertices() {  /// \todo refactor! (see grid)
      return math::ct::ipow(2u,nd);
    }

    using Vertex = NumA<nd>;
    using Vertices = std::array<Vertex,no_vertices()>;

    Ind cell_id() const { return cellId_; }
    Vertices operator()() const { return vertices_; }

    Ind cellId_;
    Vertices vertices_;
  };

  using CellVerticesRange
  = boost::transformed_range<std::function<CellVertices (NodeIdx)>, const  AnyRange<NodeIdx>>;

  /// \brief Vertices of cell at \p x_cell with cell length \p cellLength sorted
  /// in counter clockwise direction
  auto cell_vertices_coords(const Num cellLength, const NumA<nd> x_cell) const
  -> typename CellVertices::Vertices
  {
    TRACE_IN((cellLength)(x_cell));
    enum class VTK { Pixel, Polyline };
    const VTK element_type = VTK::Pixel;
    const auto cellHalfLength = 0.5 * cellLength;
    typename CellVertices::Vertices relativePosition;
    for(auto v : boost::counting_range(0U,no_edge_vertices())) {
      if(element_type == VTK::Pixel) {
        relativePosition[v] = x_cell + cellHalfLength * child_rel_pos(v).template cast<Num>();
      } else if(element_type == VTK::Polyline) {
        TERMINATE("polyline support deprecated? (TODO: decide!)");
        // old implementation, needs update:
        // relativePosition[v][d] = vertex_pos(v,d) * relLength + cellCoordinates[d];
      }
    }
    TRACE_OUT();
    return relativePosition;
  }

  /// \brief Vertices of cell at node \p nIdx sorted in counter clockwise
  /// direction
  CellVertices compute_cell_vertices(const NodeIdx nIdx) const {
    TRACE_IN((nIdx));

    const NumA<nd> x_cell = cell_coordinates(nIdx);
    const Num l_cell = cell_length(nIdx);

    return { nIdx(), cell_vertices_coords(l_cell,x_cell) };
  }

  /// \brief Total #of _leaf_ cell vertices
  /// \warning Returns the #of vertices for the leaf nodes
  /// \todo Should return the #of vertices for a grid!
  Ind no_cell_vertices() const { return nodes().no_leaf_nodes() * no_edge_vertices(); };

  /// \brief Lazy range of all _leaf_ cell vertices
  /// \warning Returns the range for the leaf nodes
  /// \todo Should return the range for a grid!
  auto cell_vertices() const -> CellVerticesRange const {
    std::function<CellVertices(const NodeIdx)> f = [&](const NodeIdx nIdx){
      return compute_cell_vertices(nIdx);
    };
    return AnyRange<Ind>(nodes().leaf_nodes()) | boost::adaptors::transformed(f);
  }

  ///@}

  /// \name Operations on mesh topology (i.e. grid connectivity)
  ///@{

  /// \brief Refines node \p nIdx isotropically.
  ///
  /// \complexity O(1) - constant time ?!? (I guess average constant time due to
  /// compression but should really check this out)
  void refine_node(const NodeIdx nIdx) {
    TRACE_IN((nIdx));
    nodes_.refine_node(nIdx);
    TRACE_OUT();
  }

  ///@}

  /// \name Operations on boundary conditions
  /// (some of this stuff is still hacky/experimental)
  ///@{

  /// \brief Range of boundaries
  inline const std::vector<Boundary>& boundaries() const { return boundaries_; }

  /// \biref Range of boundaries belonging to the solver \p solverIdx
  inline auto boundaries(const SolverIdx solverIdx) const
  -> decltype(boundaries() | boost::adaptors::filtered(std::function<bool(Boundary)>())) {
    std::function<bool(Boundary)> valid
        = [=](const Boundary& b){ return b.is_valid(solverIdx); };
    return boundaries() | boost::adaptors::filtered(valid);
  }

  /// \brief Appends boundary \p b to the grid boundaries
  void append_boundary(Boundary&& b) { boundaries_.push_back(b); }

  /// \brief Boundary cell information
  ///
  /// For each boundary in a boundary cell, the bcId as well as the solverIds
  /// that share that boundary are provided.
  struct BoundaryCell {
    BoundaryCell(const NodeIdx nodeIdx, const std::vector<SInd> boundaries)
        : nodeIdx_(nodeIdx), boundaries_(boundaries) {}
    NodeIdx node_idx() const { return nodeIdx_; }
    Ind size() const { return boundaries_.size(); }
    const std::vector<SInd>& boundaries() const { return boundaries_; }
   private:
    const NodeIdx nodeIdx_;
    const std::vector<SInd> boundaries_;
  };
  using BoundaryCells = std::vector<BoundaryCell>;

  static constexpr SInd invalid_solver(){ return invalid<SInd>(); }

  /// \brief Finds all boundary cells of \p solver
  ///
  /// Returns the subset of all cells of \p solver, that either
  /// - have no neighbor in at least one direction (grid boundary cell),
  /// or
  /// - have no neighbor belonging to solver in at least one direction
  /// (boundary cells between \p solver and other solvers).
  ///
  /// For each boundary cell, a set of boundaries is returned. For each
  /// boundary in the set, the solver ids of the corresponding solvers
  /// is returned.
  ///
  /// \warning EXPERIMENTAL!
  /// \warning Doesn't work across refinement levels yet
  ///
  BoundaryCells boundary_cells(const SolverIdx solverIdx) {
    BoundaryCells boundaryCells;

    for(auto nIdx : nodes().nodes(solverIdx)) {
      auto allCellBoundaryIds = is_cut_by_boundaries(nIdx);
      if( allCellBoundaryIds.empty() ) { continue; } // not a boundary cell
      // is a boundary cell, but maybe not from solver
      std::vector<SInd> boundaryIds; // cell might be cut by n-boundaries!
      for(auto bndryId : allCellBoundaryIds) {
        if(boundaries()[bndryId].is_valid(solverIdx)) {
          boundaryIds.push_back(bndryId);
        }
      }
      if(!boundaryIds.empty()) { // boundary cell from solver
        boundaryCells.push_back({nIdx,std::move(boundaryIds)});
      }
      // cell is not a boundary cell from solver
      // but it is a boundary cell from another solver
      // e.g this can happen over a volume coupling
    }
    return boundaryCells;
  }

  struct VolumeCoupledCell {
    Ind nodeIdx;
    std::vector<SolverIdx> solverIds;
  };

  /// \brief Finds all volume coupled cellsof \p solver
  ///
  /// Returns the subset of all cells of \p, that are also
  /// in the domain of other solvers.
  ///
  /// For each of these cells, the solver ids of all the solvers for whom they
  /// belong except the solver \p solver is returned.
  ///
  /// \warning EXPERIMENTAL
  /// \warning Doesn't work across refinement levels yet!
  template<class Solver>
  std::vector<VolumeCoupledCell> volume_coupled_cells(Solver& solver) {
    const auto solverIdx = solver.solver_idx();
    std::vector<VolumeCoupledCell> volumeCoupledCells;
    for(auto nIdx : solver.node_ids()) {
      std::vector<SInd> solverIds;
      for(auto otherSolver : nodes().grid_ids(nIdx)) {
        if(otherSolver == solverIdx) { continue; }
        solverIds.push_back(otherSolver);
      }
      if(!solverIds.empty()) {
        volumeCoupledCells.push_back(VolumeCoupledCell{solverIdx,solverIds});
      }
    }
    return volumeCoupledCells;
  }

  /// \brief EXPERIMENTAL
  Num level_set(const NodeIdx nIdx) const {
    const NumA<nd> x = cell_coordinates(nIdx);
    return level_set(x);
  }

  /// \brief EXPERIMENTAL
  Num level_set(const NumA<nd> x) const {
    Num tmp = std::numeric_limits<Num>::max();
    for(const auto& ls : boundaries_) {
      tmp = std::min(tmp,ls.signed_distance(x));
    }
    return tmp;
  }

  /// \brief EXPERIMENTAL
  template<class SignedDistance>
  bool is_cut_by(const NodeIdx nIdx, SignedDistance&& signed_distance) const {
    return is_cut_by(compute_cell_vertices(nIdx), std::forward<SignedDistance>(signed_distance));
  }

  /// \brief EXPERIMENTAL
  template<class SignedDistance>
  bool is_cut_by(const CellVertices cellVertices,
                 SignedDistance&& signed_distance) const {
    // compute level set values at vertices
    NumA<no_edge_vertices()> lsvEdges;
    int i = 0;
    for(auto vertex : cellVertices()) {
      NumA<nd> v;
      int j = 0;
      for(auto d : vertex) {
        v(j++) = d;
      }
      lsvEdges(i++) = signed_distance(v);
    }
    // if sign is not same for all then is cut!
    if((lsvEdges.array() > 0).all() || (lsvEdges.array() < 0).all()) {
      return false;
    } else {
      return true;
    }
  }

  bool is_cut_by(const CellVertices cellVertices,
                 const SInd boundaryIdx) const { // todo BoundaryIdx
    return is_cut_by(cellVertices,[&](const NumA<nd> x){
        return boundaries()[boundaryIdx].signed_distance(x);
    });
  }

  /// \brief EXPERIMENTAL
  bool is_cut_by_levelset(const NodeIdx nIdx) {
    return is_cut_by(nIdx,[&](const NumA<nd> x){ return level_set(x); });
  }

  /// \brief EXPERIMENTAL
  std::vector<SInd> is_cut_by_boundaries(const NodeIdx nIdx) {
    std::vector<SInd> result;
    SInd i = 0;
    for(auto b : boundaries()) {
      if(is_cut_by(nIdx,[&](const NumA<nd> x){ return b.signed_distance(x); })) {
        result.push_back(i);
      }
      ++i;
    }
    return result;
  }

  ///@}

  /// \brief
  /// \todo unused / deprecate?
  inline bool is_ready() const { return ready_; }

  /// \todo take grid by const reference instead of pointer (requires changes to io::Vtk!)
  friend void write_domain(const std::string fName, const This* const grid) {

    io::Vtk<nd,io::format::binary> out(grid, fName, io::precision::standard());

    out << io::stream("nodeIds",1,[](const Ind nIdx, const SInd){ return nIdx; });

    out << io::stream("nghbrIds",grid->nodes().no_samelvl_nghbr_pos(),
                      [&](const Ind nIdx, const SInd pos){
                        return grid->nodes().find_samelvl_nghbr(NodeIdx{nIdx},pos)();
                      });
    out << io::stream("solver",grid->nodes().solver_capacity(),
                      [&](const Ind nIdx, const SInd pos){
                        return grid->nodes().has_solver(NodeIdx{nIdx},SolverIdx{pos}) ? pos : invalid<SInd>();
                      });
    for(const auto& i : grid->boundaries()) {
      out << io::stream(i.name(),1,[&](const Ind nIdx, const SInd) {
          return i.signed_distance(grid->cell_coordinates(NodeIdx{nIdx}));
        });
    }
  }

  //////////////////////////////////////////////////////////////////////////////
 private:
  //////////////////////////////////////////////////////////////////////////////

  /// Contains the grid properties
  io::Properties properties_;

  /// Grid node connectivity graph
  Connectivity nodes_;

  /// Grid root cell information
  const grid::RootCell<nd> rootCell_;

  /// Domain boundaries
  std::vector<Boundary> boundaries_;

  /// Grid generator
  std::function<void(This&)> meshGeneration_;

  /// Is the grid ready to use ?
  bool ready_;

  /// \name Spatial information: implementation details
  ///@{

  /// \brief Relative position (-1 or +1) of child located at \p childPos w.r.t
  /// its parent along the \p d -axis for 2D.
  static constexpr SInt child_rel_pos_(const SInd childPos, const SInd d, D2) {
    return grid::child_rel_pos_2d_arr[childPos * nd + d];
  }

  static constexpr SIntA<nd> child_rel_pos_(const SInd childPos, D2) {
    return { child_rel_pos_(childPos,0,dim()), child_rel_pos_(childPos,1,dim()) };
  }

  /// \brief Relative position (-1 or +1) of child located at \p childPos w.r.t
  /// its parent along the \p d -axis for 3D.
  static constexpr SInt child_rel_pos_(const SInd childPos, const SInd d, D3) {
    return grid::child_rel_pos_3d_arr[childPos * nd + d];
  }

  static constexpr SIntA<nd> child_rel_pos_(const SInd childPos, D3) {
    return {child_rel_pos_(childPos,0,D3()),
            child_rel_pos_(childPos,1,D3()),
            child_rel_pos_(childPos,2,D3())};
  }

  /// \brief Relative position of neighbor in position \p nghbrPos w.r.t. cell
  /// along the \p -axis for 2D.
  static constexpr SInt nghbr_rel_pos_(const SInd nghbrPos, const SInd d, D2) {
    return grid::nghbr_rel_pos_2d_arr[nghbrPos * nd + d];
  }

  static constexpr SIntA<nd> nghbr_rel_pos_(const SInd nghbrPos, D2) {
    return {nghbr_rel_pos_(nghbrPos,0,D2()), nghbr_rel_pos_(nghbrPos,1,D2()) };
  }

  /// \brief Relative position of neighbor in position \p nghbrPos w.r.t. cell
  /// along the \p -axis for 3D.
  static constexpr SInt nghbr_rel_pos_(const SInt nghbrPos, const SInt d, D3) {
    return grid::nghbr_rel_pos_3d_arr[nghbrPos * nd + d];
  }

  static constexpr SIntA<nd> nghbr_rel_pos_(const SInd nghbrPos, D3) {
    return {nghbr_rel_pos_(nghbrPos,0,D3()),
            nghbr_rel_pos_(nghbrPos,1,D3()),
            nghbr_rel_pos_(nghbrPos,2,D3())};
  }

  /// \brief Relative position of vertex in position \p vertexId (sorted in
  /// counter clock-wise order) w.r.t. to cell along the \p -axis for 2D.
  static constexpr SInt vertex_pos_(const SInt vertexPos, const SInt d, D2) {
    return grid::vertex_pos_2d_arr[vertexPos * nd + d];
  }

  /// \brief Relative position of vertex in position \p vertexId (sorted in
  /// counter clock-wise order) w.r.t. to cell along the \p -axis for 3D.
  static constexpr SInt vertex_pos_(const SInt vertexPos, const SInt d, D3) {
    return grid::vertex_pos_3d_arr[vertexPos * nd + d];
  }

  ///@}
};

/// \brief Backwards compatibility
template<SInd nd> using Grid = CartesianHSP<nd>;

////////////////////////////////////////////////////////////////////////////////
#undef ENABLE_DBG_
////////////////////////////////////////////////////////////////////////////////
#endif
