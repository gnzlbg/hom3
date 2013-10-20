#ifndef HOM3_CONTAINER_HIERARCHICAL_IMPLEMENTATION_HPP_
#define HOM3_CONTAINER_HIERARCHICAL_IMPLEMENTATION_HPP_
////////////////////////////////////////////////////////////////////////////////
/// \file \brief Defines the hierarchical container class.
////////////////////////////////////////////////////////////////////////////////
/// Options:
#include <string>
////////////////////////////////////////////////////////////////////////////////
#define ENABLE_DBG_ 0
#include "misc/dbg.hpp"
////////////////////////////////////////////////////////////////////////////////
/// File macros:
////////////////////////////////////////////////////////////////////////////////
/// \name Assertion macros (for very common assertions)
///@{
/// \brief Asserts that a node must exist
#define assert_valid(nId) \
  ASSERT(is_valid((nId)), "Call with non-existent node.")
/// \brief Asserts that a node must be active
#define assert_active(nId) \
  ASSERT(active_((nId)), "Node " << (nId) << " is not active")
/// \brief Asserts that a child position must be in range [0,no_childs())
#define assert_child_position(pos)                              \
  ASSERT((pos) < no_child_positions(),                          \
         "Child position " << (pos) << " is not in range [0, "  \
         << no_child_positions() << ")!")
/// \brief Asserts that a child position must be in range [0,no_childs())
#define assert_neighbor_position(pos)                              \
  ASSERT((pos) < no_samelvl_neighbor_positions(),                  \
         "Neighbor position " << (pos) << " is not in range [0, "  \
         << no_samelvl_neighbor_positions() << ")!")
/// \brief Asserts that node is a leaf node
#define assert_leaf(nId) \
  ASSERT(is_leaf((nId)), "Node " << (nId) << " is not a leaf")
///@}
////////////////////////////////////////////////////////////////////////////////

namespace hom3 { namespace container {

/// \brief Hierarchical container utilities
namespace hierarchical {

namespace detail {

/// Stencil of relative sibling positions
static const constexpr std::array<std::array<SInd, 6>, 8>
rel_sibling_positions_arr = {{
  //        0               1                2                3                 4               5              << nghbrs
  {{ invalid<SInd>(),               1, invalid<SInd>(),               2, invalid<SInd>(),               4 }},  // child 0
  {{               0, invalid<SInd>(), invalid<SInd>(),               3, invalid<SInd>(),               5 }},  // child 1
  {{ invalid<SInd>(),               3,               0, invalid<SInd>(), invalid<SInd>(),               6 }},  // child 2
  {{               2, invalid<SInd>(),               1, invalid<SInd>(), invalid<SInd>(),               7 }},  // child 3
  {{ invalid<SInd>(),               5, invalid<SInd>(),               6,               0, invalid<SInd>() }},  // child 4
  {{               4, invalid<SInd>(), invalid<SInd>(),               7,               1, invalid<SInd>() }},  // child 5
  {{ invalid<SInd>(),               7,               4, invalid<SInd>(),               2, invalid<SInd>() }},  // child 6
  {{               6, invalid<SInd>(),               5, invalid<SInd>(),               3, invalid<SInd>() }}   // child 7
}};

/// Stencil of opposite neighbor positions
static const constexpr std::array<SInd, 6>
opposite_neighbor_positions_arr = {{ 1, 0, 3, 2, 5, 4 }};

/// Stencil of neighbor directions
static const constexpr std::array<SInd, 6>
neighbor_directions = {{ 0, 0, 1, 1, 2, 2 }};

}  // namespace detail

/// \brief Direction of neighbor located at \p position
inline constexpr SInd neighbor_direction(const SInd position) noexcept
{ return detail::neighbor_directions[position]; }

struct neg_dir_t {}; static const constexpr neg_dir_t neg_dir {};
struct pos_dir_t {}; static const constexpr pos_dir_t pos_dir {};

/// \brief Position of neighbor at direction \p (dir, neg_dir_t)
inline constexpr SInd neighbor_position(const SInd dir, neg_dir_t) noexcept
{ return 2 * dir; }
/// \brief Position of neighbor at direction \p (dir, pos_dir_t)
inline constexpr SInd neighbor_position(const SInd dir, pos_dir_t) noexcept
{ return 2 * dir + 1; }

/// \brief Is neighbor at \p position in the negative side of its direction?
inline constexpr bool is_neighbor_at(const SInd position, neg_dir_t) noexcept
{ return position % 2 == 0; }

/// \brief Is neighbor at \p position in the positive side of its direction?
inline constexpr bool is_neighbor_at(const SInd position, pos_dir_t) noexcept
{ return position % 2 != 0; }

static inline constexpr
SInd no_samelvl_neighbor_positions(const SInd nd) noexcept
{ return 2 * nd; }

/// \brief Hierarchical Node Container: stores a hierarchical "nd"-space
/// partition that is able to store "ng" grids.
///
/// Notation:
///
/// Idx  -> Index
/// Ids  -> Indices
/// nIdx -> nodeIdx
/// pIdx -> parentIdx
/// cIdx -> childIdx
/// pos  -> position
/// cPos -> child position (in parent)
/// nPos -> neighbor position (relative to node)
///
/// Child positions: child(nId,pos)
/// <pre>
///                  / z(3)
///                /
///            6 | 7
///            -----
///            4 | 5
///          /
///   ^ y (1)
///   |   /
/// 2 | 3
/// --|-----> x (0)
/// 0 | 1
/// </pre>
///
/// Face/Neighbor positions: nghbr(nId,pos)
/// <pre>
///            ^ 5
///           /
///       3 /
///     _____
///  0 | nId | 1
///    ------
///   /   2
/// 4
/// </pre>
///
template<SInd nd_> struct Implementation {
  static const SInd nd = nd_;  ///< \f$\#\f$ of space dimensions
  using This = Implementation<nd>;

  //////////////////////////////////////////////////////////////////////////////
  /// \name Encoded spatial information
  ///
  /// Some information is encoded in the position of data within arrays
  /// The following functions give the \f$\#\f$ of positions
  ///@{

  /// \brief \f$\#\f$ of child positions
  static inline constexpr SInd no_child_positions()
  { return math::ct::ipow(2u, nd); }
  /// \brief \f$\#\f$ of same level neighbor positions
  static inline constexpr SInd no_samelvl_neighbor_positions()
  { return container::hierarchical::no_samelvl_neighbor_positions(nd); }
  /// \brief Maximum number of allowed levels
  static inline constexpr SInd max_no_levels() { return 20; }
  ///@}
  //////////////////////////////////////////////////////////////////////////////

  /// \brief Constructs a grid container that can maximally hold \p maxNoNodes
  /// nodes and \p maxNoSolvers grids.
  ///
  /// \param [in] maxNoGridNodes   maximum number of nodes that the container
  ///                              can hold.
  /// \param [in] maxNoGridSolvers maximum number of solver grids that the
  ///                              container can store.
  Implementation(io::Properties input)
      : noActiveNodes_{0}
      , maxNoNodes_{io::read<Ind>(input, "maxNoGridNodes")}
      , lowerFreeNodeBound_{0}
      , parentIds_{this, "parents"}
      , childrenIds_{this, "childs"}
      , isFree_{this, "isFree"}
      , node2cells_{maxNoNodes_, io::read<SInd>(input, "maxNoGridSolvers")} {
    TRACE_IN_();
    initialize_root_node_();
    TRACE_OUT();
  }
  Implementation() = delete;
  std::string name() const noexcept { return "hierarchical"; }

  //////////////////////////////////////////////////////////////////////////////
  /// \name Status operations
  ///@{

  /// \brief \f$\#\f$ of nodes that the container can allocate
  inline Ind  capacity()      const noexcept { return maxNoNodes_;         }
  /// \brief \f$\#\f$ of nodes in use
  inline Ind  size()          const noexcept { return noActiveNodes_;      }
  /// \brief nId of the first node in the tree
  inline NodeIdx node_begin() const noexcept { return NodeIdx{0};          }
  /// \brief nId of the last node _in use_ in the tree
  inline NodeIdx node_end()   const noexcept { return NodeIdx{size()};     }
  /// \brief nId of the last node in the tree
  inline NodeIdx node_last()  const noexcept { return NodeIdx{capacity()}; }
  /// \brief \f$\#\f$ of nodes that are not free
  inline Ind  no_nodes()      const noexcept { return noNodes_;            }
  /// \brief Is the container empty?
  inline bool empty()         const noexcept { return noActiveNodes_ == 0; }

  ///@}
  //////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////////
  /// \name Node-to-solver-cell mapping
  ///@{

  /// \brief Cell idx of node \p nodeIdx within \p solverIdx grid
  inline CellIdx cell_idx
  (const NodeIdx nodeIdx, const SolverIdx solverIdx) const noexcept {
    assert_valid(nodeIdx); assert_active(nodeIdx);
    return node2cells_(nodeIdx(), solverIdx());
  }
  /// \brief Reference to the cell idx of node \p nodeIdx within \p solverIdx
  /// grid
  inline CellIdx& cell_idx
  (const NodeIdx nodeIdx, const SolverIdx solverIdx)      noexcept {
    assert_valid(nodeIdx); assert_active(nodeIdx);
    return node2cells_(nodeIdx(), solverIdx());
  }
  /// \brief Is the node \p nodeIdx part of \p solverIdx 's grid ?
  inline bool has_solver
  (const NodeIdx nodeIdx, const SolverIdx solverIdx) const noexcept {
    assert_valid(nodeIdx); assert_active(nodeIdx);
    return is_valid(cell_idx(nodeIdx, solverIdx));
  }
  /// \brief \f$\#\f$ of solver grids that the container can hold
  inline SInd solver_capacity() const noexcept { return node2cells_.cols(); }

  ///@}
  //////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////////
  /// \name Node operations
  ///@{

  /// \brief Is node \p nIdx an active node?
  /// (active means "in use", i.e. not free)
  inline bool active_(const NodeIdx nIdx) const noexcept
  { assert_valid(nIdx); return !isFree_(nIdx()); }
  /// \brief Is the node \p nIdx a root node?
  /// (right now there is only a single root node)
  inline bool is_root(const NodeIdx nIdx) const noexcept
  { assert_valid(nIdx); return !is_valid(parent(nIdx)); }
  /// \brief Is node \p nIdx a leaf node?
  inline bool is_leaf(const NodeIdx nIdx) const noexcept
  { assert_valid(nIdx); assert_active(nIdx); return no_childs(nIdx) == 0; }

  /// \brief \f$\#\f$ of leaf nodes in the container
  inline Ind no_leaf_nodes() const noexcept
  { return boost::distance(leaf_nodes()); }

  /// \brief Parent of node \p nIdx
  inline NodeIdx parent(const NodeIdx nIdx) const noexcept
  { assert_valid(nIdx); return parentIds_(nIdx()); }

  /// \brief Child of \p nIdx at position \p pos
  inline NodeIdx child(const NodeIdx nIdx, const SInd pos) const noexcept {
    assert_valid(nIdx); assert_child_position(pos);
    return childrenIds_(nIdx()) + NodeIdx{pos};
  }
  /// \brief \f$\#\f$ of children of \p nIdx
  inline SInd no_childs(const NodeIdx nIdx) const noexcept
  { assert_valid(nIdx); return boost::distance(childs(nIdx)); }
  /// \brief Does node \p nIdx have children?
  inline bool has_children(const NodeIdx nIdx) const noexcept
  { assert_valid(nIdx); return !is_leaf(nIdx); }
  /// \brief Does node \p nIdx have _all_ children? I.e.
  /// is \f$\#\f$ of children of \p nIdx == \f$\#\f$ of children positions?
  inline bool has_all_children(const NodeIdx nIdx) const noexcept
  { assert_valid(nIdx); return no_childs(nIdx) == no_child_positions(); }
  /// \brief Position in parent of node \p nIdx
  SInd position_in_parent(const NodeIdx nIdx) const noexcept {
    TRACE_IN((nIdx));
    assert_valid(nIdx); assert_active(nIdx);
    const auto pIdx = parent(nIdx);
    assert_valid(pIdx); assert_active(pIdx);
    const SInd res = boost::find_if(child_positions(), [&](const SInd cPos) {
        return child(pIdx, cPos) == nIdx; }) - child_positions().begin();
    TRACE_OUT();
    return res;
  }

  /// \brief \f$\#\f$ of same level neighbors of node \p nIdx
  inline SInd no_samelvl_neighbors(const NodeIdx nIdx) const noexcept
  { assert_valid(nIdx); return boost::distance(samelvl_neighbors(nIdx)); }

  /// \brief Opposite neighbor position. I.e. the neighbor with
  /// position \p pos of a given nIdx, has nIdx as neighbor in its opposite
  /// position.
  ///
  /// \warning returns iSInd if out of bounds. Bound check
  /// not performed here.
  static constexpr SInd opposite_neighbor_position(const SInd pos) noexcept {
    assert_neighbor_position(pos);
    return detail::opposite_neighbor_positions_arr[pos];
  }

  /// \brief Position in parent of a sibling in direction "nghbrPos"
  /// from a sibling in position "childPos"
  ///
  /// \param [in] childPos position in parent of a node id
  /// \param [in] nghbrPos neighbor position with respect to nodePos
  ///
  /// \returns position in common parent of a sibling located at nghbrPos w.r.t.
  /// nodePos if it exist
  static constexpr SInd rel_sibling_position
  (const SInd childPos, const SInd nghbrPos) noexcept {
    assert_child_position(childPos); assert_neighbor_position(nghbrPos);
    return detail::rel_sibling_positions_arr[childPos][nghbrPos];
  }
  ///@}
  //////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////////
  /// \name Modifying algorithms
  ///@{

  /// \brief Refines a leaf node isotropically, i.e. into 4 nodes in 2D or into
  /// 8 nodes in 3D, and returns the node idx of the first children.
  ///@{
  NodeIdx refine_node(const NodeIdx nIdx) noexcept {
    TRACE_IN((nIdx)); using namespace algorithm;
    assert_active(nIdx); assert_leaf(nIdx);
    ASSERT(level(nIdx) + 1 < max_no_levels(),
           "can't refine > max_no_levels()!");

    const auto firstChildIdx = is_compact() ? node_end() : free_spot_();
    const auto lastChildIdx = firstChildIdx + NodeIdx{no_child_positions()};

    no_nodes_() += no_child_positions();
    lowerFreeNodeBound_ += NodeIdx{no_child_positions()};
    //size_() = size() + no_child_positions();
    if(lastChildIdx >= node_end()) {
      auto diff = math::absdiff(lastChildIdx, node_end());
      while(diff != NodeIdx{0}) {
        ++size_();
        --diff;
      }
    }

    ASSERT([&]() {
      return all_of(Range<NodeIdx>(firstChildIdx, lastChildIdx),
                    [&](const NodeIdx cIdx) { return !is_free(cIdx); }); }(),
      "All future childrens must be free!");

    ASSERT([&]() {
      return all_of(Range<NodeIdx>(firstChildIdx, lastChildIdx),
                    [&](const NodeIdx cIdx) { return !is_reseted(cIdx);} ); }(),
      "All future childrens must be reseted!");

    child_(nIdx) = firstChildIdx;
    for (const auto& childIdx : childs(nIdx)) {
      parent_(childIdx) = nIdx;
      isFree_(childIdx()) = false;
    }

    ASSERT([&]() {  // Assert that the children exist and are fine
        for (const auto& cIdx : childs(nIdx)) {
          assert_valid(cIdx); assert_active(cIdx); assert_leaf(cIdx);
        }
        return true;
      }(), "post-condition on newly created children");
    TRACE_OUT();
    return firstChildIdx;
  }

  /// \brief Inserts a node into parent \p pIdx at position \p pos.
  ///
  /// \complexity O(1) if the tree is_compact()
  /// \complexity average O(1) if the tree is not Compact
  // NodeIdx insert_node(const NodeIdx pIdx, const SInd pos) {
  //   TERMINATE("unimplemented");
  //   // TRACE_IN((pIdx)(pos));
  //   // assert_active(pIdx);
  //   // ASSERT(!is_valid(child(pIdx)), "position must be free!");
  //   // ASSERT(level(pIdx) + 1 < max_no_levels(),
  //   //        "can't refine > max_no_levels()!");
  //   // const auto nIdx = is_compact() ? node_end() : free_spot_();
  //   // ++no_nodes_();
  //   // ++lowerFreeNodeBound_;
  //   // if(nIdx == node_end()) {
  //   //   ++size_();
  //   // }
  //   // reset_node_(nIdx);
  //   // child_(pIdx,pos) = nIdx;
  //   // parent_(nIdx) = pIdx;
  //   // isFree_(nIdx()) = false;
  //   // TRACE_OUT();
  //   // return nIdx;
  // }

  /// \brief Removes node \p nIdx from its parent.
  ///
  /// \complexity: O(1)
  // NodeIdx remove_node(const NodeIdx /*nIdx*/) {
  //   TERMINATE("unimplemented");
  //   // TRACE_IN((nIdx));
  //   // assert_active(nIdx);
  //   // const auto pIdx = parent(nIdx);
  //   // ASSERT(is_valid(pIdx), "Node has no parent!");
  //   // const auto pos = position_in_parent(nIdx);
  //   // ASSERT(is_valid(pos), "Node is not child of parent!");
  //   // child_(pIdx,pos) = invalid<Ind>();
  //   // reset_node_(nIdx);
  //   // --no_nodes_();
  //   // if(nIdx == size()) {
  //   //   --size();
  //   // }
  //   // isReady_ = false;
  //   // lowerFreeNodeBound_ = std::min(lowerFreeNodeBound_,nIdx);
  //   // TRACE_OUT();
  //   // return pIdx;
  // }

  ///@}

  ///@}
  //////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////////
  /// \name Non-modifying algorithms
  ///@{

  /// \brief Computes the level of the node \p nIdx
  /// \complexity
  ///   - O(\f$\#\f$l): linear in the number of levels,
  ///   - O(log(\f$\#\f$size)): logarithmic in the number of cells
  inline SInd level(NodeIdx nIdx) const noexcept {
    TRACE_IN((nIdx));
    assert_active(nIdx);
    SInd level_ = 0;
    nIdx = parent(nIdx);
    while (is_valid(nIdx)) {  // traverse the tree up to the root cell
      assert_active(nIdx);
      nIdx = parent(nIdx);
      ++level_;  // counting the levels
    }
    ASSERT(level_ < max_no_levels(), "error, level " + std::to_string(level_)
           + " exceeds max level " + std::to_string(max_no_levels()) + "!");
    TRACE_OUT();
    return level_;
  }

  /// \brief Computes the same level neighbor of the node \p nIdx located at
  /// its relative position \p pos.
  ///
  /// \param [in] nIdx  node id of the cell whose neighbor is to be determined.
  /// \param [in] pos  position of the required neighbor w.r.t. the node nIdx.
  ///
  /// \returns [Ind] node idx of the neighbor located at pos w.r.t. the node
  /// nIdx. If no neighbor is found it returns invalid<Ind>().
  ///
  /// \complexity unknown (worst case is probably O(\f$\#\f$l) = 2 * \f$\#\f$l,
  /// need to traverse the tree up to the root and back down again).
  /// \todo Improve complexity and profile!
  ///
  /// \requires a balanced tree that satisfies 2:1 rule.
  NodeIdx find_samelvl_neighbor
  (const NodeIdx nIdx, const SInd nghbrPos) const noexcept {
    TRACE_IN((nIdx)(nghbrPos)); using namespace algorithm;
    DBG("start samelvl_neighbor | nIdx: ", nIdx, " | nghbrPos ", nghbrPos);
    assert_valid(nIdx); assert_active(nIdx); assert_neighbor_position(nghbrPos);

    /// The root cell has no neighbors:
    if (is_root(nIdx)) {
      DBG("nIdx ", nIdx, " is a root cell | no nghbr found in pos: ", nghbrPos);
      TRACE_OUT();
      return invalid<NodeIdx>();
    }

    // Vector of traversed positions
    memory::stack::arena<SInd, max_no_levels()> stackMemory;
    auto traversedPositions = memory::stack::make<std::vector>(stackMemory);

    /// Travel up the tree until a node i is found, s.t. the previously
    /// traversed node is a child of i in the position that is "opposite" to pos
    /// (opposite in the sense of opposite_neighbor_position(nghbrPos), see that
    /// function call for details). The traversed positions are saved.
    auto commonParentFound = false;
    auto currentNode = nIdx;
    DBG("starting up-traversal at node: ", nIdx);
    while (!commonParentFound) {
      DBGV((currentNode));
      const auto pIdx = parent(currentNode);
      const auto posInParent = position_in_parent(currentNode);
      DBGV((parent(currentNode))(posInParent)
           (rel_sibling_position(posInParent, nghbrPos)));

      /// If the parent has a sibling in direction "nghbrPos" we have found
      /// a common parent and are done
      if (is_valid(rel_sibling_position(posInParent, nghbrPos))) {
        DBG("common parent found -> up-traversal finished");
        commonParentFound = true;
        break;
      }

      if (is_root(pIdx)) {
        /// If the parent is the root and there is no sibling in direction
        /// "nbghrPos", cell has no neighbor and we are done
        DBG("parent is root, and no sibling in direction found!",
            "-> nghbr not found for nIdx: ", nIdx, ", nghbrPos: ", nghbrPos);
        TRACE_OUT();
        return invalid<NodeIdx>();
      } else {
        /// Keep on going up the tree
        traversedPositions.emplace_back(posInParent);
        currentNode = pIdx;
      }
    }
    ASSERT(commonParentFound, "No common parent has been found!");

    /// This is the opposite node at the common parent:
    currentNode = child(parent(currentNode),
                        rel_sibling_position(position_in_parent(currentNode),
                                             nghbrPos));

    DBG("starting down traversal at node: ", currentNode);
    /// Traverse the tree back from the opposite node in opposite order to find
    /// the neighbor:
    for (const auto& traversedPosition : traversedPositions | reversed) {
      const auto nextChildPosition
          = rel_sibling_position(traversedPosition,
                                 opposite_neighbor_position(nghbrPos));
      const auto nextIdx = child(currentNode, nextChildPosition);
      DBGV((currentNode)(traversedPosition)(nextChildPosition)(nextIdx));
      if (is_valid(nextIdx)) {
        currentNode = nextIdx;
      } else {
        DBG("nextIdx: ", nextIdx, " does not exist!",
            " -> nghbr not found for nIdx ", nIdx, ", nghbrPos: " , nghbrPos);
        TRACE_OUT();
        return invalid<NodeIdx>();
      }
    }

    DBG("nghbr found for nIdx: ", nIdx, ", nghbrPos: ", nghbrPos, "!",
        " -> nghbrIdx: " , currentNode);

    assert_valid(currentNode);  // Current node must be a valid neighbor
    ASSERT(level(nIdx) == level(currentNode),
           "Neighbors are not at the same level!");
    TRACE_OUT();
    return currentNode;
  }

  /// \brief Returns the same level neighbors of the node nIdx
  ///
  /// \param [in] nIdx node id of the cell whose neighbors are to be computed.
  /// \returns [IndA<no_samelvl_nghbr_positions()>] Array containing the same
  /// level neighbors of node \p nIdx.
  NodeIdxA<no_samelvl_neighbor_positions()>
  find_samelvl_neighbors(const NodeIdx nIdx) const noexcept {
    TRACE_IN((nIdx));
    assert_valid(nIdx); assert_active(nIdx);
    NodeIdxA<no_samelvl_neighbor_positions()> nodeNghbrs;
    for (const auto& pos : neighbor_positions()) {
      nodeNghbrs[pos] = find_samelvl_neighbor(nIdx, pos);
    }
    TRACE_OUT();
    return nodeNghbrs;
  }

  /// \brief Returns in which neighbor position \p nghbrIdx is w.r.t. the node
  /// \p nIdx (returns invalid position if not found).
  SInd which_neighbor(const NodeIdx nIdx,
                      const NodeIdx nghbrNodeIdx) const noexcept {
    SInd pos = 0;
    for (const auto& nghbr : all_samelvl_neighbors(nIdx)) {
      if (nghbr == nghbrNodeIdx) {
        return pos;
      }
      ++pos;
    }
    pos = invalid<SInd>();  // not found, returns invalid position

    // Consistency check:
    ASSERT([&]() {
        using std::to_string;
        const auto conjugate_pos = which_neighbor(nghbrNodeIdx, nIdx);
        if (pos != invalid<SInd>()) {
          // pos == valid, conjugate which_neighbor should return opposite pos
          ASSERT(conjugate_pos == opposite_neighbor_position(pos),
                 "nIdx = " + to_string(nIdx()) +
                 " doesn't have nghbrNodeIdx = " + nghbrNodeIdx
                 + " in the opposite neighbor position "
                 + opposite_neighbor_position(pos) + " but in the position "
                 + conjugate_pos + "!");
        } else {
          // pos == invalid, conjugate which_neighbor should return invalid
          ASSERT(conjugate_pos == pos,
                 "nIdx = " + to_string(nIdx)
                 + " doesn't have nghbrNodeIdx = " + nghbrNodeIdx
                 + " as neighbor, but nghbrNodeIdx has nIdx as neighbor in pos "
                 + conjugate_pos + "!");
        }
      }(), "consistency check");
    return pos;
  }


  ///@}
  //////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////////
  /// \name Node range filters
  ///@{

  /// \brief Returns [RangeFilter] of existing node ids
  inline auto valid() const noexcept -> RangeFilter<NodeIdx> {
    return {[&](const NodeIdx nIdx) { return is_valid(nIdx); }};
  }
  /// \brief Returns [RangeFilter] of active node ids
  inline auto active() const noexcept  -> RangeFilter<NodeIdx> {
    return {[&](const NodeIdx nIdx) { return active_(nIdx); }};
  }
  /// \brief Returns [RangeFilter] of leaf node ids:
  inline auto leafs() const noexcept -> RangeFilter<NodeIdx> {
    return {[&](const NodeIdx nIdx) {return is_leaf(nIdx);}};
  }
  /// \brief Returns [RangeFilter] of node ids at level "l"
  inline auto atLevel(const SInd l) const noexcept -> RangeFilter<NodeIdx> {
    return {[&, l](const NodeIdx nIdx) { return level(nIdx) == l; }};
  }
  /// \brief Maps child positions to child indices
  inline auto child_position_to_idx(const NodeIdx nIdx) const noexcept
  -> RangeTransformer<SInd, NodeIdx> {
    return {[&, nIdx](const SInd pos) {
        assert_valid(nIdx); assert_child_position(pos);
        const auto firstChildIdx = child_(nIdx);
        // if cell is leaf, child_(nIdx) is invalid == numeric_limits<Ind>::max
        // -> adding stuff to it makes it either wrap around or is UB!
        return is_valid(firstChildIdx) ? firstChildIdx + NodeIdx{pos}
                                       : invalid<NodeIdx>(); }};
  }
  /// \brief Maps neighbor positions to neighbor indices
  inline auto neighbor_position_to_idx(const NodeIdx nIdx) const noexcept
  -> RangeTransformer<SInd, NodeIdx> {
    return {[&, nIdx](const SInd pos) {
              assert_valid(nIdx); assert_neighbor_position(pos);
              return find_samelvl_neighbor(nIdx, pos); }};
  }

  ///@}
  //////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////////
  /// \name Node ranges
  ///@{

  /// \brief Returns [IndRange] of _all_ node position Idxs. That is,
  /// including inactive nodes.
  inline auto all_nodes() const noexcept -> Range<NodeIdx>
  { return {node_begin(), node_end()}; }
  /// \brief Returns [FilteredRange] of all active node Idxs.
  inline auto nodes() const RETURNS(all_nodes() | active());
  /// \brief Returns [FilteredRange] of all active leaf node Idxs.
  ///
  /// \todo Remove return type.
  inline auto leaf_nodes() const RETURNS(nodes() | leafs());
  /// \brief Returns [IndRange] of all child positions.
  inline auto child_positions() const noexcept -> Range<SInd>
  { return {SInd(0), no_child_positions()}; }
  /// \brief Returns [IndRange] of _all_ child Idxs of node \p nIdx.
  /// That is, including inactive childs.
  inline auto all_childs(const NodeIdx nIdx) const
  RETURNS(child_positions() | child_position_to_idx(nIdx));
  /// \brief Returns [FilteredRange] of existing child Idxs of node \p nIdx.
  ///
  /// \todo Remove return type.
  inline auto childs(const NodeIdx nIdx) const
  RETURNS(all_childs(nIdx) | valid() | active());

  /// \brief Returns [IndRange] of all neighbor positions.
  inline auto neighbor_positions() const noexcept -> Range<SInd>
  { return {SInd{0}, no_samelvl_neighbor_positions()}; }
  /// \brief Returns [Range] of _all_ neighbor Ids of node nIdx. That is,
  /// including inactive neighbors.
  ///
  /// \todo Remove return type.
  inline auto all_samelvl_neighbors(const NodeIdx nIdx, strict) const
  RETURNS(find_samelvl_neighbors(nIdx));
  inline auto all_samelvl_neighbors(const NodeIdx nIdx, lazy) const
  RETURNS(neighbor_positions() | neighbor_position_to_idx(nIdx));
  template<class T = lazy>
  inline auto all_samelvl_neighbors(const NodeIdx nIdx, T t = T()) const
  RETURNS(all_samelvl_neighbors(nIdx, t));

  /// \brief Returns [FilteredRange] of existing neighbor Ids of node \p nIdx.
  ///
  /// \todo Remove return type.
  inline auto samelvl_neighbors(const NodeIdx nIdx) const
  RETURNS(all_samelvl_neighbors(nIdx) | valid());

  ///@}
  //////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////////
  /// \name Solver grids cell ranges
  ///@{

  /// \brief Filters nodes that belong to the \p solverIdx 's grid
  inline auto in_solver(const SolverIdx solverIdx) const noexcept
  -> RangeFilter<NodeIdx> {
    return {[&, solverIdx](const NodeIdx nIdx) {
              return has_solver(nIdx, solverIdx); }};
  }

  /// \brief Transforms a range of solver ids into a range of cells ids within
  /// each solver's grid for the node \p nIdx
  inline auto solver_position_to_cell(const NodeIdx nIdx) const noexcept
  -> RangeTransformer<SolverIdx, CellIdx> {
    return {[&, nIdx](const SolverIdx solverIdx) {
              return cell_idx(nIdx, solverIdx); }};
  }

  /// \brief Transforms a range of nodes into a range of cell ids belonging to
  /// the \p solverIdx 's grid
  inline auto node_to_cell(const SolverIdx solverIdx) const noexcept
  -> RangeTransformer<NodeIdx, CellIdx> {
    return {[&, solverIdx](const NodeIdx nIdx) {
        return is_valid(nIdx) ? cell_idx(nIdx, solverIdx)
                              : invalid<CellIdx>(); }};
  }

  /// \brief Returns range of _all_ solver ids.
  inline auto solver_ids() const noexcept -> Range<SolverIdx>
  { return {SolverIdx{0}, SolverIdx{solver_capacity()}}; }

  /// \brief Returns range of _all_ cell ids of node \p nIdx,
  /// including invalid ones!
  ///
  /// \todo Remove return type.
  inline auto all_cell_ids(const NodeIdx nIdx) const
  RETURNS(solver_ids() | solver_position_to_cell(nIdx));

  /// \brief Return range of all nodes belonging to the \p solverIdx 's grid
  ///
  /// \todo Remove return type.
  inline auto nodes(const SolverIdx solverIdx) const
  RETURNS(nodes() | in_solver(solverIdx));

  /// \brief Return range of all cells belonging to the \p solverIdx 's grid
  ///
  /// \todo Remove return type.
  inline auto cells(const SolverIdx solverIdx) const
  RETURNS(nodes(solverIdx) | node_to_cell(solverIdx));

  /// \brief Range of all same level neighbors of \p nIdx belonging to \p sIdx
  inline auto all_samelvl_neighbors(const NodeIdx nIdx, SolverIdx sIdx) const
  RETURNS(all_samelvl_neighbors(nIdx) | node_to_cell(sIdx));

  ///@}
  //////////////////////////////////////////////////////////////////////////////


 private:
  //////////////////////////////////////////////////////////////////////////////
  /// \name Container Status Implementation details
  ///@{

  ///     ________________________________________________
  ///    |   |   |   |   |   |   |   |   |   |   |   |   | noActiveNodes = 6
  ///    | 0 | 1 | x | 2 | x | 3 | 4 | 5 | x | x | x | x | noNodes = 8
  ///    |___|___|___|___|___|___|___|___|___|___|___|___| maxNoNodes = 12
  /// first^   empty^    lastNodeInUse^        capacity^

  Ind noActiveNodes_;  ///< \f$\#\f$ of active nodes
  Ind maxNoNodes_;     ///< max \f$\#\f$ of nodes
  Ind noNodes_;        ///< \f$\#\f$ of nodes (all node idx > noNodes_ are
                       ///< not in use)

  NodeIdx lowerFreeNodeBound_;  ///< Smallest free node id

  //////////////////////////////////////////////////////////////////////////////
  /// \name Node Data
  ///@{

  /// Alias for data containers:
  template<template <SInd> class V_, SInd nd__ = 1>
  using M = container::Matrix
    <Implementation, container::matrix::tag::Cell, V_, Ind, SInd, nd__>;


  M<NodeIdxM> parentIds_;    ///< Parent node ids for each node
  M<NodeIdxM> childrenIds_;  ///< Child node ids for each node
  M<BoolMatrix> isFree_;     ///< Indicates if a node is free or in use

  ///@}
  //////////////////////////////////////////////////////////////////////////////

  /// \brief Is the \f$\#\f$ of nodes equalto the \f$\#\f$ of active nodes?
  bool is_compact() const noexcept { return noActiveNodes_ == noNodes_; };

  /// Mapping from node indices to solver cell indices (one for each solver
  /// grid)
  EigenDynRowMajor<CellIdx> node2cells_;

  /// \brief Smallest node id that is not in use. If all
  /// node ids are in use, returns size (i.e. one past the end).
  inline NodeIdx free_spot_() const noexcept {
    const auto spot
      = boost::find_if(Range<NodeIdx>(lowerFreeNodeBound_, node_end()),
                       [&](const NodeIdx i) { return isFree_(i()); }) - std::begin(all_nodes());
    return NodeIdx{spot};
  }

  /// \brief Writable reference to noActiveNodes_
  inline Ind& size_() noexcept { return noActiveNodes_; }
  /// \brief Writable reference no noNodes_
  inline Ind& no_nodes_() noexcept { return noNodes_; }

  /// \brief Is a node position in memory free?
  /// (is that position not holding a node?)
  ///
  /// \todo !is_valid(parent(nIdx)) && !is_valid(child(nIdx));
  inline bool is_free(const NodeIdx nIdx) const noexcept {
    assert_valid(nIdx);
    return isFree_(nIdx());
  }

  ///@}
  //////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////////
  /// \name Write access to container variables
  ///@{

  /// \brief Writable reference to parent id of node id \p nIdx
  inline NodeIdx& parent_(const NodeIdx nIdx) noexcept
  { assert_valid(nIdx); return parentIds_(nIdx()); }
  /// \brief Writable reference to child id of \p nIdx at position \p pos
  inline NodeIdx& child_(const NodeIdx nIdx) noexcept
  { assert_valid(nIdx); return childrenIds_(nIdx()); }
  inline NodeIdx child_(const NodeIdx nIdx) const noexcept
  { assert_valid(nIdx); return childrenIds_(nIdx()); }

  ///@}
  //////////////////////////////////////////////////////////////////////////////


  //////////////////////////////////////////////////////////////////////////////
  /// \name Initialization routines
  ///@{

  /// \brief Creates an initial root node
  ///
  /// \warning The container must be empty!
  void initialize_root_node_() {
    TRACE_IN_();
    ASSERT(empty(), "Container is not empty!");
    ++size_(); ++no_nodes_();
    isFree_(node_begin()()) = false;  // activate node before reseting it
    reset_node_(node_begin());
    isFree_(node_begin()()) = false;
    TRACE_OUT();
  }

  /// \brief Resets data of all nodes
  void reset_nodes_(const NodeIdx first, const NodeIdx last) {
    TRACE_IN_();
    for (const auto& nIdx : Range<NodeIdx>(first, last)) {
      reset_node_(nIdx);
    }
    TRACE_OUT();
  }

  /// \brief Resets all the ids in \p nIdx to invalid<Ind>() and marks node as
  /// free.
  ///
  /// NOTE: need to keep in sync with is_reseted !
  void reset_node_(const NodeIdx nIdx) {
    TRACE_IN((nIdx));
    parent_(nIdx) = invalid<NodeIdx>();
    child_(nIdx) = invalid<NodeIdx>();
    for (const auto& pos : solver_ids()) {
      cell_idx(nIdx, pos) = invalid<CellIdx>();
    }
    isFree_(nIdx()) = true;
    TRACE_OUT();
  }

  ///@}
  //////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////////
  /// \name Debugging methods
  ///@{

  /// \brief Is a node reseted?
  bool is_reseted(const NodeIdx nIdx) {
    using namespace algorithm;
    TRACE_IN((nIdx));
    TRACE_OUT();
    return !is_valid(parent(nIdx))
        && !is_valid(child_(nIdx))
        && all_of(all_cell_ids(nIdx),
                  [&](const CellIdx cIdx) { return !is_valid(cIdx); })
        && is_free(nIdx);
  }

  ///@}
  //////////////////////////////////////////////////////////////////////////////
};

}  // namespace hierarchical

template<SInd nd> using Hierarchical = hierarchical::Implementation<nd>;

}  // namespace container

////////////////////////////////////////////////////////////////////////////////
}  // namespace hom3
////////////////////////////////////////////////////////////////////////////////
#undef assert_valid
#undef assert_active
#undef assert_child_position
#undef ENABLE_DBG_
////////////////////////////////////////////////////////////////////////////////
#endif
