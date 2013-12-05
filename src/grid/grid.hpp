#ifndef HOM3_GRID_HPP_
#define HOM3_GRID_HPP_
////////////////////////////////////////////////////////////////////////////////
/// \file \brief Contains the grid data-structure
///
/// \todo Replace warnings by compile-time/run-time assertions
////////////////////////////////////////////////////////////////////////////////
/// Includes:
#include <limits>
#include <algorithm>
#include <vector>
#include "globals.hpp"
#include "containers/hierarchical.hpp"
#include "generation.hpp"
#include "boundary.hpp"
#include "io/output.hpp"
#include <iomanip>
#include <Eigen/Geometry>
/// Options:
#define ENABLE_DBG_ 0
#include "misc/dbg.hpp"
////////////////////////////////////////////////////////////////////////////////

namespace hom3 { namespace grid {

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

static constexpr std::array<std::tuple<SInd,Num>,4*2> vertex_vertex_connectivity_2d_arr {{
    // x          y        z
    {1, +1.}, {2, +1.}, // pos: 0
    {0, -1.}, {3, +1.}, // pos: 1
    {3, +1.}, {0, -1.}, // pos: 2
    {2, -1.}, {1, -1.}  // pos: 3
}};

static constexpr std::array<std::tuple<SInd,Num>,8*3> vertex_vertex_connectivity_3d_arr {{
    // x          y        z
    {1, +1.}, {2, +1.}, {4, +1.},  // pos: 0
    {0, -1.}, {3, +1.}, {5, +1.},  // pos: 1
    {3, +1.}, {0, -1.}, {6, +1.},  // pos: 2
    {2, -1.}, {1, -1.}, {7, +1.},  // pos: 3
    {5, +1.}, {6, +1.}, {0, -1.},  // pos: 4
    {4, -1.}, {7, +1.}, {1, -1.},  // pos: 5
    {7, +1.}, {4, -1.}, {2, -1.},  // pos: 6
    {6, -1.}, {5, -1.}, {1, -1.}   // pos: 7
}};

///@}


template<int nd> struct RootCell {
  RootCell(const RootCell<nd>&) = default;
  RootCell(const NumA<nd>& min, const NumA<nd>& max) : length(math::eps) {
    ASSERT([&](){
      NumA<nd> lengths = max - min;
      Num previous_length = lengths(0);
      for (SInd d = 1; d != nd; ++d) {
        ASSERT(math::approx(lengths(d), previous_length),
          "Error length mismatch between dimensions d = "
          + std::to_string(d - 1) + " (x_min = " + std::to_string(min(d - 1))
          + ", x_max = " + std::to_string(max(d - 1))
          + ", length = " + std::to_string(lengths(d - 1)) + "),"
          + " and d = " + std::to_string(d) + " (x_min = " + std::to_string(min(d))
          + ", x_max = " + std::to_string(max(d))
          + ", length = " + std::to_string(lengths(d)) + ")."
          + " Root cell is not square shaped!");
        }
        return true;
      }(), "Root cell must be square shaped!");

    for (SInd d = 0; d != nd; ++d) {
      const Num lengthD = max(d) - min(d); // move up! length( (max-min).max() )
      ASSERT(lengthD > math::eps, "Negative length not allowed!");
      length = std::max(length, lengthD);
      coordinates(d) = min(d) + 0.5 * lengthD; // move up! coords( min + 0.5 * length)
    }
  }

  static const int nDim = nd;
  Num length;
  NumA<nd> coordinates;
  NumA<nd> x_min() const noexcept { return coordinates.array() - 0.5 * length; }
  NumA<nd> x_max() const noexcept { return coordinates.array() + 0.5 * length; }
};

/// \brief Hierarchical Cartesian Grid data-structure
/// - couples the spatial information (e.g. coordinates) with the connectivity
/// graph (e.g. neighbor/parent-child relationships)
/// - manages the boundaries/geometry of each sub domain
///
/// \todo rename to CartesianHSP
template <SInd nd_> struct CartesianHSP : container::Hierarchical<nd_> {

  static const SInd nd = nd_;
  using This = CartesianHSP<nd>;
  using Boundary = boundary::Interface<nd>;
  using Boundaries = std::vector<Boundary>;
  using Dim = dim<nd>;
  using D2 = dim<2>;
  using D3 = dim<3>;

  /// \name Connectivity requirements (for internal use)
  ///@{
  using Connectivity = typename container::Hierarchical<nd>;
  using Connectivity::nodes;
  using Connectivity::leaf_nodes;
  using Connectivity::no_leaf_nodes;
  using Connectivity::level;
  using Connectivity::parent;
  using Connectivity::position_in_parent;
  ///@}

  /// \name Construction
  ///@{

  /// \brief Construct a grid from a set of input properties
  CartesianHSP(io::Properties input)
      : container::Hierarchical<nd>(input)
      , properties_(input)
      , rootCell_(io::read<RootCell<nd>>(properties_,"rootCell"))
      , ready_(false)
  { TRACE_IN_(); TRACE_OUT(); }

  /// \brief Constructs and initializes a grid from a set of input properties
  CartesianHSP(io::Properties input, initialize_t) : CartesianHSP(input) {
    TRACE_IN_();
    initialize();
    TRACE_OUT();
  }

  /// \brief Initializes the grid, i.e. reads the boundaries and the mesh
  /// generator from the properties if they haven't been set yet and generates
  /// the mesh
  void initialize() {
    //ASSERT(!boundaries().empty(), "Grid has no boundaries!"); // -> warning!
    if(!meshGeneration_) {
      read_mesh_generator();
      generate_mesh();
    }
    ready_ = true;
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
  inline auto dimensions() const -> Range<SInd> { return {SInd{0}, nd}; }
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
    return child_rel_pos_(childPos,d,Dim());
  }

  /// \brief Relative position of child localted at \p childPos w.r.t. its parent
  ///
  /// \returns vector of relative positions
  static constexpr SIntA<nd> child_rel_pos(const SInd childPos) {
    return child_rel_pos_(childPos,Dim());
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
    return nghbr_rel_pos_(nghbrPos,d,Dim());
  }

  /// \brief Relative position of neighbor localted at \p nghbrPos w.r.t. its parent
  ///
  /// \returns vector of relative positions
  static constexpr SIntA<nd> nghbr_rel_pos(const SInd nghbrPos) {
    return nghbr_rel_pos_(nghbrPos,Dim());
  }

  /// \brief Relative position of vertex in position \p vertexId (sorted in
  /// counter clock-wise order) w.r.t. to cell along the \p -axis.
  ///
  /// That is:
  ///              o________________________o
  ///            /|   pos: 6   |   pos: 7  /|
  ///           / | (-1,+1,+1) | (+1,+1,+1) |
  ///          /  |____________|____________|
  ///         /   |   pos: 4   |   pos: 5   |
  ///        /    | (-1,-1,+1) | (+1,-1,+1) |
  ///       /     o____________|____________o
  ///      /     /                   /     /
  ///     o_____/__________________o/     /
  ///    |   pos: 2   |   pos: 3   |     /    y (1) ^
  ///    | (-1,+1,-1) | (+1,+1,-1) |    /           |     ^ z (2)
  ///    |____________|____________|   /            |    /
  ///    |   pos: 0   |   pos: 1   |  /             |  /
  ///    | (-1,-1,-1) | (+1,-1,-1) | /              |/
  ///    o____________|____________o/               o-------> x (0)
  ///
  /// Note:
  ///  - the positions correspond to the vertices.
  ///  - the stencil is the same one as that for the relative child position!
  ///  - the ordering agrees with VTK's Pixel (8) and Voxel(11) elements!
  static constexpr SIntA<nd> vertex_rel_pos(const SInd vertexIdx) {
    return child_rel_pos_(vertexIdx, Dim());
  }

  /// \brief Length of cells at \p level
  Num cell_length_at_level(const Ind level) const {
    return rootCell_.length / std::pow(2,level);
  }

  /// \brief Length of cell at node \p nIdx
  Num cell_length(const NodeIdx nIdx) const {
    return cell_length_at_level(level(nIdx));
  }

  /// \brief Centroid coordinates of cell at node \p nIdx
  ///
  /// Note: the coordinates are computed from the root cell's coordinates.
  /// \complexity O(L) - linear in the number of levels.
  /// \complexity O(logN) - logarithmic in the number of cells.
  NumA<nd> cell_coordinates(const NodeIdx nIdx) const {
    TRACE_IN((nIdx));
    const Ind level_ = level(nIdx);
    if( level_ != 0 ) {
      const auto pIdx      = parent(nIdx);
      const auto posInParent  = position_in_parent(nIdx);
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

  /// \brief Centroid coordinates of neighbor at \p nghbrPos of node \p nIdx
  ///
  /// Works even if the neighbor doesn't exist.
  NumA<nd> neighbor_coordinates(const NodeIdx nIdx, const SInd nghbrPos) const {
    TRACE_IN((nIdx)(nghbrPos));
    const auto x_node = cell_coordinates(nIdx);
    const auto length = cell_length(nIdx);

    NumA<nd> x_nghbr;
    for (auto d : dimensions()) {
      x_nghbr(d) = x_node(d) + nghbr_rel_pos(nghbrPos,d) * length;
    }
    TRACE_OUT();
    return x_nghbr;
  }

  /// \todo rename cell_id to node_id
  struct CellVertices {
    static constexpr SInd no_vertices() {  /// \todo refactor! (see grid)
      return math::ct::ipow(2u,nd);
    }

    using Vertex = NumA<nd>;
    using Vertices = std::array<Vertex,no_vertices()>;

    Ind cell_idx() const { return cellIdx_; }
    Vertices operator()() const { return vertices_; }

    Ind cellIdx_;
    Vertices vertices_;
  };

  using CellVerticesRange
  = boost::transformed_range<std::function<CellVertices(NodeIdx)>,
                             const AnyRange<NodeIdx>>;

  auto vertex_positions() const
  RETURNS(boost::counting_range(SInd{0},no_edge_vertices()));

  /// \brief Vertices of cell at \p x_cell with cell length \p cellLength
  /// Note: the order is same as the order of a cells children
  auto cell_vertices_coords(const Num cellLength, const NumA<nd> x_cell) const
  -> typename CellVertices::Vertices {
    TRACE_IN((cellLength)(x_cell));
    const auto cellHalfLength = 0.5 * cellLength;
    typename CellVertices::Vertices relativePosition;
    for(auto v : vertex_positions()) {
      relativePosition[v] = x_cell + cellHalfLength
                            * vertex_rel_pos(v).template cast<Num>();
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
  Ind no_cell_vertices() const
  { return no_leaf_nodes() * no_edge_vertices(); };

  /// \brief Lazy range of all _leaf_ cell vertices
  /// \warning Returns the range for the leaf nodes
  /// \todo Should return the range for a grid!
  auto cell_vertices() const -> CellVerticesRange {
    std::function<CellVertices(const NodeIdx)> f = [&](const NodeIdx nIdx){
      return compute_cell_vertices(nIdx);
    };
    return AnyRange<NodeIdx>(leaf_nodes()) | boost::adaptors::transformed(f);
  }

  ///@}

  /// \name Operations on boundary conditions
  ///@{

  /// \brief Appends boundary \p b to the grid boundaries
  void append_boundary(Boundary b) noexcept {
    boundaries_.emplace_back(std::move(b));
  }

  /// \brief Range of _all_ boundaries
  inline const Boundaries& boundaries() const { return boundaries_; }

  /// \biref Range of \p solverIdx 's boundaries
  inline auto boundaries(const SolverIdx solverIdx) const
  -> decltype(boundaries()
              | boost::adaptors::filtered(std::function<bool(Boundary)>())) {
    std::function<bool(Boundary)> valid
        = [=](const Boundary& b){ return b.is_valid(solverIdx); };
    return boundaries() | boost::adaptors::filtered(valid);
  }


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

    for (auto nIdx : nodes(solverIdx)) {
      auto allCellBoundaryIds = is_cut_by_boundaries(nIdx);
      if( allCellBoundaryIds.empty() ) { continue; } // not a boundary cell
      // is a boundary cell, but maybe not from solver
      std::vector<SInd> boundaryIds; // cell might be cut by n-boundaries!
      for (auto bndryId : allCellBoundaryIds) {
        if(boundaries()[bndryId].solver_idx() == solverIdx) {
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

  /// \brief Filters those nodes cut by the \p boundary .
  template<class Boundary>
  inline auto cut_by_boundary(Boundary&& boundary) -> RangeFilter<NodeIdx> {
    return {[&, boundary](const NodeIdx nIdx) {
        return is_cut_by(nIdx, [&, boundary](const NumA<nd> x) {
            return boundary.signed_distance(x); }
          ); }
    };
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
      for(auto otherSolver : grid_ids(nIdx)) {
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
  template<class SignedDistance>
  bool is_cut_by(const NodeIdx nIdx, SignedDistance&& signed_distance) const {
    auto lsvEdges = signed_distance_at_vertices(nIdx, std::forward<SignedDistance>(signed_distance));
    // if sign is not same for all then is cut!
    if((lsvEdges.array() > 0).all() || (lsvEdges.array() < 0).all()) {
      return false;
    } else {
      return true;
    }
  }

  template<class SignedDistance>
  NumA<no_edge_vertices()> signed_distance_at_vertices(const NodeIdx nIdx,
                                                       SignedDistance&& signed_distance) const {
    const auto vertices = compute_cell_vertices(nIdx)();
    NumA<no_edge_vertices()> sd_at_vertices;
    for (auto pos : vertex_positions()) {
      sd_at_vertices(pos) = signed_distance(vertices[pos]);
    }
    return sd_at_vertices;
  }

  template<class Vertices, class LSVs, class CutPoints>
  void debug_cut_points_at_node(const NodeIdx nIdx, const Vertices& vertices,
                                const LSVs& lsv,
                                const CutPoints& cutPoints) const {

    auto print_point = [&](const std::string label, auto vertex) {
      std::cerr << label;
      for (auto d : dimensions()) { std::cerr << vertex(d) << " "; }
      std::cerr << "\n";
    };
    auto cp = [](auto c) { return std::get<0>(c); };
    auto vB_p = [](auto c) { return std::get<1>(c); };
    auto vB = [&](auto c) { return vertices()[vB_p(c)]; };
    auto vE_p = [](auto c) { return std::get<2>(c); };
    auto vE = [&](auto c) { return vertices()[vE_p(c)]; };

    std::cerr << std::setprecision(16);
    std::cerr << "Wrong #of cut points: " << cutPoints.size()
              <<  " in node: " << nIdx << ".\n";

    print_point("Cell center coordinates are:\n  x_c: ", cell_coordinates(nIdx));

    std::cerr << "\n";
    std::cerr << "Cell length is: " << cell_length(nIdx) << ".\n";
    std::cerr << "Vertices are:\n";
    for (auto v : vertex_positions()) {
      print_point("  #" + std::to_string(v) + ", x: ", vertices()[v]);
      std::cerr << " lsv: " << lsv[v] << "\n";
    }
    std::cerr << "Cut-points are:\n";
    SInd i = 0;
    for (auto c : cutPoints) {
      std::cerr << "  #" << i << " vB_p: " << vB_p(c) << " vE_p: " << vE_p(c) << "\n";
      print_point("    vB_x: ", vB(c));
      print_point("    vE_x: ", vE(c));
      print_point("    CP_x: ", cp(c));
      ++i;
    }
  }

  struct CutSurface {
    NumA<nd> centroid;
    std::array<NumA<nd>, nd> cut_points;
    NumA<nd> normal;
  };

  template<class SignedDistance>
  CutSurface cut_surface(const NodeIdx nIdx, SignedDistance&& signed_distance) const {
    using namespace geometry::algorithm;
    using std::get;
    ASSERT(is_valid(nIdx), "invalid node!");
    ASSERT(is_cut_by(nIdx, signed_distance), "node must be cut by sd!");

    const auto vertices = compute_cell_vertices(nIdx);
    std::array<Num, no_edge_vertices()> vertices_lsv;
    for (auto v : vertex_positions()) {
      vertices_lsv[v] = signed_distance(vertices()[v]);
    }
    memory::stack::arena<std::tuple<NumA<nd>,SInd,SInd>, no_edge_vertices()+2> stackMemory;
    auto cutPoints = memory::stack::make<std::vector>(stackMemory);

    for (auto v : vertex_positions()) {
      for (auto d : dimensions()) {
        SInd ov; Num s;
        std::tie(ov, s) = vertex_vertex_connectivity(v, d);
        ASSERT(!math::approx(vertices_lsv[v], 0.) && !math::approx(vertices_lsv[ov], 0.),
               "lsv at vertex is zero!");
        if (math::signum(vertices_lsv[v]) != math::signum(vertices_lsv[ov])) {
          const auto cp = cut_point_along_direction<nd>
                          (vertices()[v], vertices()[ov],
                           vertices_lsv[v], vertices_lsv[ov], d);
          cutPoints.emplace_back(std::make_tuple(std::move(cp), v, ov));
        }
      }
    }

    remove_duplicate_points<nd>(cutPoints,
      [&](auto a, auto b) { return PointLT<nd>()(get<0>(a), get<0>(b)); },
      [&](auto a, auto b) { return PointEQ()(get<0>(a), get<0>(b)); });

    ASSERT([&](){
        if((nd == 2 && cutPoints.size() == 2)
           || (nd == 3 && (cutPoints.size() == 3 || cutPoints.size() == 4))) {
          return true;
        }
        debug_cut_points_at_node(nIdx, vertices, vertices_lsv, cutPoints);
        return false;
      }(), "wrong number of cutpoints!");

    auto centroid = surface_centroid<nd>()([&](const SInd p){ return std::get<0>(cutPoints[p]); });
    auto normal = surface_normal<nd>()([&](const SInd p){ return std::get<0>(cutPoints[p]); });

    auto scaled_normal = Eigen::Scaling(cell_length(nIdx) * 0.25) * normal;
    auto translated_centroid = Eigen::Translation<Num, nd>(scaled_normal) * centroid;
    //DBGV_ON((normal)(centroid)(scaled_normal)(translated_centroid)(signed_distance(translated_centroid)));
    if (signed_distance(translated_centroid) < 0.) {
      normal = Eigen::Rotation2D<Num>(math::pi) * normal;
    }

    CutSurface cs;
    cs.centroid = centroid;
    cs.cut_points[0] = std::get<0>(cutPoints[0]);
    cs.cut_points[1] = std::get<0>(cutPoints[1]);
    cs.normal = normal;
    return cs;
  }

  template<class SignedDistance>
  NumA<nd> cut_surface_centroid(const NodeIdx nIdx, SignedDistance&& signed_distance) const {
    using namespace geometry::algorithm;
    using std::get;
    ASSERT(is_valid(nIdx), "invalid node!");
    ASSERT(is_cut_by(nIdx, signed_distance), "node must be cut by sd!");

    const auto vertices = compute_cell_vertices(nIdx);
    std::array<Num, no_edge_vertices()> vertices_lsv;
    for (auto v : vertex_positions()) {
      vertices_lsv[v] = signed_distance(vertices()[v]);
    }
    memory::stack::arena<std::tuple<NumA<nd>,SInd,SInd>, no_edge_vertices()+2> stackMemory;
    auto cutPoints = memory::stack::make<std::vector>(stackMemory);

    for (auto v : vertex_positions()) {
      for (auto d : dimensions()) {
        SInd ov; Num s;
        std::tie(ov, s) = vertex_vertex_connectivity(v, d);
        ASSERT(!math::approx(vertices_lsv[v], 0.) && !math::approx(vertices_lsv[ov], 0.),
               "lsv at vertex is zero!");
        if (math::signum(vertices_lsv[v]) != math::signum(vertices_lsv[ov])) {
          const auto cp = cut_point_along_direction<nd>
                          (vertices()[v], vertices()[ov],
                           vertices_lsv[v], vertices_lsv[ov], d);
          cutPoints.emplace_back(std::make_tuple(std::move(cp), v, ov));
        }
      }
    }

    remove_duplicate_points<nd>(cutPoints,
      [&](auto a, auto b) { return PointLT<nd>()(get<0>(a), get<0>(b)); },
      [&](auto a, auto b) { return PointEQ()(get<0>(a), get<0>(b)); });

    ASSERT([&](){
        if((nd == 2 && cutPoints.size() == 2)
           || (nd == 3 && (cutPoints.size() == 3 || cutPoints.size() == 4))) {
          return true;
        }
        debug_cut_points_at_node(nIdx, vertices, vertices_lsv, cutPoints);
        return false;
      }(), "wrong number of cutpoints!");

    return surface_centroid<nd>()([&](const SInd p){ return std::get<0>(cutPoints[p]); });
  }


  bool is_cut_by(const CellVertices cellVertices,
                 const SInd boundaryIdx) const { // todo BoundaryIdx
    return is_cut_by(cellVertices,[&](const NumA<nd> x){
        return boundaries()[boundaryIdx].signed_distance(x);
    });
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
  bool is_cut_by_levelset(const NodeIdx nIdx) {
    return is_cut_by(nIdx, [&](const NumA<nd> x){ return level_set(x); });
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

  /// \brief Root cell
  RootCell<nd> root_cell() const noexcept { return rootCell_; }

  friend void write_domain(const This& grid) {
    write_domain(std::to_string(nd) + "D_Grid", grid);
  }

  /// \todo take grid by const reference instead of pointer (requires changes to
  /// io::Vtk!)
  friend void write_domain(const String fName, const This& grid) {
    io::Vtk<nd, io::format::ascii> out(&grid, fName, io::precision::standard());

    out << io::stream("nodeIds", 1, [](const Ind nIdx, const SInd) {
      return nIdx;
    });

    out << io::stream("nghbrIds", grid.no_samelvl_neighbor_positions(),
      [&](const Ind nIdx, const SInd pos) {
        return grid.find_samelvl_neighbor(NodeIdx{nIdx}, pos)();
    });
    out << io::stream("x_center", grid.no_dimensions(),
      [&](const Ind nIdx, const SInd pos) {
        return grid.cell_coordinates(NodeIdx{nIdx})(pos);
    });
    out << io::stream("solver", grid.solver_capacity(),
      [&](const Ind nIdx, const SInd pos) {
        return grid.has_solver(NodeIdx{nIdx}, SolverIdx{pos}) ?
               pos : invalid<SInd>();
    });
    for(const auto& i : grid.boundaries()) {
      out << io::stream(i.name(), 1, [&](const Ind nIdx, const SInd) {
        return i.signed_distance(grid.cell_coordinates(NodeIdx{nIdx}));
      });
    }
  }

  //////////////////////////////////////////////////////////////////////////////
 private:
  //////////////////////////////////////////////////////////////////////////////

  /// Contains the grid properties
  io::Properties properties_;

  /// Grid root cell information
  const RootCell<nd> rootCell_;

  /// Domain boundaries
  Boundaries boundaries_;

  /// Grid generator
  std::function<void(This&)> meshGeneration_;

  /// Is the grid ready to use ?
  bool ready_;

  /// \name Spatial information: implementation details
  ///@{

  /// \brief Relative position (-1 or +1) of child located at \p childPos w.r.t
  /// its parent along the \p d -axis for 2D.
  static constexpr SInt child_rel_pos_(const SInd childPos, const SInd d, D2) {
    return child_rel_pos_2d_arr[childPos * nd + d];
  }

  static constexpr SIntA<nd> child_rel_pos_(const SInd childPos, D2) {
    return { child_rel_pos_(childPos,0,Dim()), child_rel_pos_(childPos,1,Dim()) };
  }

  /// \brief Relative position (-1 or +1) of child located at \p childPos w.r.t
  /// its parent along the \p d -axis for 3D.
  static constexpr SInt child_rel_pos_(const SInd childPos, const SInd d, D3) {
    return child_rel_pos_3d_arr[childPos * nd + d];
  }

  static constexpr SIntA<nd> child_rel_pos_(const SInd childPos, D3) {
    return {child_rel_pos_(childPos,0,D3()),
            child_rel_pos_(childPos,1,D3()),
            child_rel_pos_(childPos,2,D3())};
  }

  /// \brief Relative position of neighbor in position \p nghbrPos w.r.t. cell
  /// along the \p -axis for 2D.
  static constexpr SInt nghbr_rel_pos_(const SInd nghbrPos, const SInd d, D2) {
    return nghbr_rel_pos_2d_arr[nghbrPos * nd + d];
  }

  static constexpr SIntA<nd> nghbr_rel_pos_(const SInd nghbrPos, D2) {
    return {nghbr_rel_pos_(nghbrPos,0,D2()), nghbr_rel_pos_(nghbrPos,1,D2()) };
  }

  /// \brief Relative position of neighbor in position \p nghbrPos w.r.t. cell
  /// along the \p -axis for 3D.
  static constexpr SInt nghbr_rel_pos_(const SInt nghbrPos, const SInt d, D3) {
    return nghbr_rel_pos_3d_arr[nghbrPos * nd + d];
  }

  static constexpr SIntA<nd> nghbr_rel_pos_(const SInd nghbrPos, D3) {
    return {nghbr_rel_pos_(nghbrPos,0,D3()),
            nghbr_rel_pos_(nghbrPos,1,D3()),
            nghbr_rel_pos_(nghbrPos,2,D3())};
  }

  static constexpr std::tuple<SInd, Num>
  vertex_vertex_connectivity_(const SInd i, D2)
  { return vertex_vertex_connectivity_2d_arr[i]; }

  static constexpr std::tuple<SInd, Num>
  vertex_vertex_connectivity_(const SInd i, D3)
  { return vertex_vertex_connectivity_3d_arr[i]; }

  /// Vertex sharing \p edge with vertex at position \p vertexPos
  /// Returns tuple of (vertex_position, sign)
  static constexpr std::tuple<SInd, Num>
  vertex_vertex_connectivity(const SInd vertexPos, const SInt d) {
    return vertex_vertex_connectivity_(vertexPos * nd + d, Dim());
  }

  ///@}
};

/// \brief Backwards compatibility
template<SInd nd> using Grid = CartesianHSP<nd>;

////////////////////////////////////////////////////////////////////////////////
}} // hom3::grid namespace
////////////////////////////////////////////////////////////////////////////////
#undef ENABLE_DBG_
////////////////////////////////////////////////////////////////////////////////
#endif
