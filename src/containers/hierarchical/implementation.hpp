#ifndef HOM3_CONTAINER_HIERARCHICAL_IMPLEMENTATION_HPP_
#define HOM3_CONTAINER_HIERARCHICAL_IMPLEMENTATION_HPP_
////////////////////////////////////////////////////////////////////////////////
/// \file \brief Defines the hierarchical container class.
////////////////////////////////////////////////////////////////////////////////
/// Options:
#define ENABLE_DBG_ 0
#include "../../misc/dbg.hpp"
////////////////////////////////////////////////////////////////////////////////
/// File macros:
////////////////////////////////////////////////////////////////////////////////
/// \name Assertion macros (for very common assertions)
///@{
/// \brief Asserts that a node must exist
#define assert_valid(nId) \
  ASSERT(is_valid((nId)),"Call with non-existent node.")
/// \brief Asserts that a node must be active
#define assert_active(nId) \
  ASSERT(active_((nId)),"Node " << (nId) << " is not active")
/// \brief Asserts that a child position must be in range [0,no_childs())
#define assert_child_position(pos) \
  ASSERT((pos) < no_child_pos(), "Child position " << (pos) << " does not exist!")
///@}
////////////////////////////////////////////////////////////////////////////////

namespace hom3 { namespace container {

/// \brief Hierarchical container utilities
namespace hierarchical {

/// Stencil of relative sibling positions
static constexpr std::array<std::array<SInd,6>,8>  rel_sibling_position_arr {{
  //        0       1       2        3         4        5   | nghbrPos
  {{ invalid<SInd>(),               1, invalid<SInd>(),               2, invalid<SInd>(),               4 }}, // child 0
  {{               0, invalid<SInd>(), invalid<SInd>(),               3, invalid<SInd>(),               5 }}, // child 1
  {{ invalid<SInd>(),               3,               0, invalid<SInd>(), invalid<SInd>(),               6 }}, // child 2
  {{               2, invalid<SInd>(),               1, invalid<SInd>(), invalid<SInd>(),               7 }}, // child 3
  {{ invalid<SInd>(),               5, invalid<SInd>(),               6,               0, invalid<SInd>() }}, // child 4
  {{               4, invalid<SInd>(), invalid<SInd>(),               7,               1, invalid<SInd>() }}, // child 5
  {{ invalid<SInd>(),               7,               4, invalid<SInd>(),               2, invalid<SInd>() }}, // child 6
  {{               6, invalid<SInd>(),               5, invalid<SInd>(),               3, invalid<SInd>() }}  // child 7
}};

/// Stencil of opposite neighbor positions
static constexpr std::array<SInd,6> opposite_nghbr_position_arr{{ 1, 0, 3, 2, 5, 4 }};


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
template<SInd nd> struct Implementation {

  typedef Implementation<nd> This;

  //////////////////////////////////////////////////////////////////////////////
  /// \name Encoded spatial information
  ///
  /// Some information is encoded in the position of data within arrays
  /// The following functions give the \f$\#\f$ of positions
  ///@{

  /// \brief \f$\#\f$ of child positions
  static inline constexpr SInd no_child_pos() { return math::ct::ipow(2u,nd); }
  /// \brief \f$\#\f$ of same level neighbor positions
  static inline constexpr SInd no_samelvl_nghbr_pos() { return 2 * nd; }
  /// \brief Maximum number of allowed levels
  static inline constexpr SInd max_no_levels() { return 20; }
  ///@}
  //////////////////////////////////////////////////////////////////////////////

  /// \brief Constructs a grid container that can maximally hold \p maxNoNodes
  /// nodes and \p maxNoSolvers grids.
  ///
  /// \param [in] maxNoNodes maximum number of nodes that the container can hold.
  /// \param [in] maxNoSolvers maximum number of solver grids that the container can store.
  Implementation(const Ind maxNoNodes, const SInd maxNoSolvers)
      : noActiveNodes_{0}, maxNoNodes_{maxNoNodes}, lowerFreeNodeBound_{0},
    parentIds_(this,"parents"), childrenIds_(this,"childs"), isFree_(this,"isFree"),
    isReady_{false}, node2cells_{maxNoNodes,maxNoSolvers} {
    TRACE_IN((maxNoNodes)(maxNoSolvers));

    initialize_root_node_();

    TRACE_OUT();
  }
  Implementation() = delete;
  std::string name() const { return "hierarchical"; }

  //////////////////////////////////////////////////////////////////////////////
  /// \name Status operations
  ///@{

  /// \brief \f$\#\f$ of nodes that the container can allocate
  inline Ind  capacity() const { return maxNoNodes_;         }
  /// \brief \f$\#\f$ of nodes in use
  inline Ind  size()     const { return noActiveNodes_;      }
  /// \brief nId of the first node in the tree
  inline NodeIdx node_begin() const { return NodeIdx{0}; }
  /// \brief nId of the last node _in use_ in the tree
  inline NodeIdx node_end()   const { return NodeIdx{size()}; }
  /// \brief nId of the last node in the tree
  inline NodeIdx node_last()  const { return NodeIdx{capacity()}; }
  /// \brief \f$\#\f$ of nodes that are not free
  inline Ind  no_nodes() const { return noNodes_;            }
  /// \brief Is the container empty?
  inline bool empty()    const { return noActiveNodes_ == 0; }
  /// \brief \f$\#\f$ of space dimensions of the container
  inline SInd nDim()     const { return nd;                  }
  /// \brief Is the container ready for usage? (no modifications are currently
  /// taking place, all invariants hold).
  inline SInd is_ready() const { return isReady_;            }

  ///@}
  //////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////////
  /// \name Node-to-solver-cell mapping
  ///@{

  /// \brief Cell id of node \p nodeIdx within \p solverIdx grid
  inline const CellIdx& cell_id(const NodeIdx nodeIdx, const SolverIdx solverIdx) const {
    assert_valid(nodeIdx); assert_active(nodeIdx);
    return node2cells_(nodeIdx(),solverIdx());
  }
  /// \brief Reference to the cell id of node \p nodeIdx within \p solverIdx grid
  inline CellIdx& cell_id(const NodeIdx nodeIdx, const SolverIdx solverIdx) {
    assert_valid(nodeIdx); assert_active(nodeIdx);
    return node2cells_(nodeIdx(),solverIdx());
  }
  /// \brief Is the node \p nodeIdx part of \p solverIdx 's grid ?
  inline bool has_solver(const NodeIdx nodeIdx, const SolverIdx solverIdx) const {
    assert_valid(nodeIdx); assert_active(nodeIdx);
    return is_valid(cell_id(nodeIdx,solverIdx));
  }
  /// \brief \f$\#\f$ of solver grids that the container can hold
  inline SInd solver_capacity() const {
    return node2cells_.cols();
  }

  ///@}
  //////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////////
  /// \name Node operations
  ///@{

  /// \brief Is node \p nIdx an active node?
  /// (active means "in use", i.e. not free)
  inline bool active_(const NodeIdx nIdx) const {
    TRACE_IN((nIdx));
    assert_valid(nIdx);
    TRACE_OUT();
    return !isFree_(nIdx());
  }
  /// \brief Is the node \p nIdx a root node?
  /// (right now there is only a single root node)
  inline bool is_root(const NodeIdx nIdx) const {
    TRACE_IN((nIdx));
    assert_valid(nIdx);
    TRACE_OUT();
    return !is_valid(parent(nIdx));
  }
  /// \brief Is node \p nIdx a leaf node?
  inline bool is_leaf(const NodeIdx nIdx) const {
    TRACE_IN((nIdx));
    assert_valid(nIdx); assert_active(nIdx);
    TRACE_OUT();
    return no_childs(nIdx) == 0;
  }
  /// \brief \f$\#\f$ of leaf nodes in the container
  inline Ind no_leaf_nodes() const { TRACE_IN_(); TRACE_OUT();
    return boost::distance(leaf_nodes());
  }

  /// \brief Parent of node \p nIdx
  inline const NodeIdx& parent(const NodeIdx nIdx) const { assert_valid(nIdx); return parentIds_(nIdx()); }

  /// \brief Child of \p nIdx at position \p pos
  inline const NodeIdx& child(const NodeIdx nIdx, const SInd pos) const {
    TRACE_IN((nIdx)(pos));
    assert_valid(nIdx); assert_child_position(pos);
    TRACE_OUT();
    return childrenIds_(nIdx(),pos);
  }
  /// \brief \f$\#\f$ of children of \p nIdx
  inline SInd no_childs(const NodeIdx nIdx) const {
    TRACE_IN((nIdx));
    assert_valid(nIdx);
    TRACE_OUT();
    return boost::distance(childs(nIdx));
  }
  /// \brief Does node \p nIdx have children?
  inline bool has_children(const NodeIdx nIdx) const { return !is_leaf(nIdx); }
  /// \brief Does node \p nIdx have _all_ children? I.e.
  /// is \f$\#\f$ of children of \p nIdx == \f$\#\f$ of children positions?
  inline bool has_all_children(const NodeIdx nIdx) const {
    TRACE_IN((nIdx));
    assert_valid(nIdx);
    TRACE_OUT();
    return no_childs(nIdx) == no_child_pos();
  }
  /// \brief Position in parent of node \p nIdx
  SInd position_in_parent(const NodeIdx nIdx) const {
    TRACE_IN((nIdx));
    assert_valid(nIdx); assert_active(nIdx);
    const auto pIdx = parent(nIdx);
    assert_valid(pIdx); assert_active(pIdx);
    const SInd res = boost::find_if(child_pos(), [&](const SInd cPos){
        return child(pIdx,cPos) == nIdx; }) - child_pos().begin();
    TRACE_OUT();
    return res;
  }

  /// \brief \f$\#\f$ of same level neighbors of node \p nIdx
  inline SInd no_samelvl_nghbrs(const NodeIdx nIdx) const {
    TRACE_IN((nIdx));
    assert_valid(nIdx);
    TRACE_OUT();
    return boost::distance(samelvl_nghbrs(nIdx));
  }

  /// \brief Opposite nghbr position. I.e. the neighbor with
  /// position \p pos of a given nIdx, has nIdx as neighbor in its opposite
  /// position.
  ///
  /// \warning returns iSInd if out of bounds. Bound check
  /// not performed here.
  static constexpr SInd opposite_nghbr_position(const SInd pos) {
    return DBG_EXPR(pos < 6 ?)
        opposite_nghbr_position_arr[pos]
        DBG_EXPR(: invalid<SInd>());
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
  (const SInd childPos, const SInd nghbrPos) {
    return DBG_EXPR(childPos < 8 && nghbrPos < 6 ?)
           rel_sibling_position_arr[childPos][nghbrPos]
           DBG_EXPR(: invalid<SInd>());
  }
  ///@}
  //////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////////
  /// \name Modifying algorithms
  ///@{

  /// \brief Inserts a node into parent \p pIdx at position \p pos.
  ///
  /// \complexity O(1) if the tree is_compact()
  /// \complexity average O(1) if the tree is not Compact
  NodeIdx insert_node(const NodeIdx pIdx, const SInd pos) {
    TRACE_IN((pIdx)(pos));
    assert_active(pIdx);
    assert_valid(!child(pIdx,pos)); // position must be free!
    ASSERT(level(pIdx) + 1 < max_no_levels(), "can't refine > max_no_levels()!");
    const auto nIdx = is_compact() ? node_end() : free_spot_();
    ++no_nodes_();
    ++lowerFreeNodeBound_;
    if(nIdx == node_end()) {
      ++size_();
    }
    reset_node_(nIdx);
    child_(pIdx,pos) = nIdx;
    parent_(nIdx) = pIdx;
    isFree_(nIdx()) = false;
    TRACE_OUT();
    return nIdx;
  }

  /// \brief Removes node \p nIdx from its parent.
  ///
  /// \complexity: O(1)
  NodeIdx remove_node(const NodeIdx nIdx) {
    TRACE_IN((nIdx));
    assert_active(nIdx);
    const auto pIdx = parent(nIdx);
    ASSERT(is_valid(pIdx), "Node has no parent!");
    const auto pos = position_in_parent(nIdx);
    ASSERT(is_valid(pos), "Node is not child of parent!");
    child_(pIdx,pos) = invalid<Ind>();
    reset_node_(nIdx);
    --no_nodes_();
    if(nIdx == size()) {
      --size();
    }
    isReady_ = false;
    lowerFreeNodeBound_ = std::min(lowerFreeNodeBound_,nIdx);
    TRACE_OUT();
    return pIdx;
  }

  /// \brief Refines a node isotropically, i.e. into 4 (in 2D) or 8 (in 3D)
  /// children. It optionally takes a predicate \p predicate (const SInd
  /// posInParent) that is executed after each child is inserted.
  ///@{
  template<class P> NodeIdx refine_node(const NodeIdx nIdx, P&& predicate) {
    TRACE_IN((nIdx));
    for(auto pos : child_pos()) {
      insert_node(nIdx,pos); // note insert_node returns nIdx of child
      predicate(pos);        // but computing it from the nIdx is easy with pos.
    }
    isReady_ = false;
    TRACE_OUT();
    return child(nIdx,0);
  }
  NodeIdx refine_node(const NodeIdx nIdx) { return refine_node(nIdx,[](const SInd){}); }
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
  inline SInd level(NodeIdx nIdx) const {
    TRACE_IN((nIdx));
    assert_active(nIdx);
    SInd level_ = 0;
    nIdx = parent(nIdx);
    while(is_valid(nIdx)) {
      assert_active(nIdx);
      nIdx = parent(nIdx);
      ++level_;
    }
    ASSERT(level_ < max_no_levels(), "error!");
    TRACE_OUT();
    return level_;
  }

  /// \brief Computes the same level neighbor of the node \p nIdx located at
  /// its relative position \p pos.
  ///
  /// \param [in] nIdx  node id of the cell whose neighbor is to be determined.
  /// \param [in] pos  position of the required neighbor w.r.t. the node nIdx.
  ///
  /// \returns [Ind] node id of the neighbor located at pos w.r.t. the node nIdx.
  /// If no neighbor is found it returns invalid<Ind>().
  ///
  /// \complexity unknown (worst case is probably O(\f$\#\f$l) = 2 * \f$\#\f$l,
  /// need to traverse the tree up to the root and back down again).
  /// \todo Improve complexity and profile!
  ///
  /// \requires a balanced tree that satisfies 2:1 rule.
  NodeIdx find_samelvl_nghbr(const NodeIdx nIdx, const SInd nghbrPos) const {
    TRACE_IN((nIdx)(nghbrPos));
    DBG("start samelvl_nghbr | nIdx: ", nIdx, " | nghbrPos ", nghbrPos);

    /// The root cell has no neighbors:
    if(is_root(nIdx)) {
      DBG("nIdx ", nIdx, " is a root cell | no nghbr found in nghbrPos: ", nghbrPos);
      TRACE_OUT();
      return invalid<NodeIdx>();
    }

    // Create vector of reserved positions on the stack
    memory::stack::arena<SInd,max_no_levels()> stackMemory;
    auto traversedPositions = memory::stack::make<std::vector>(stackMemory);
    traversedPositions.reserve(max_no_levels() /*level(nIdx)*/);

    /// Travel up the tree until a node i is found, s.t. the previously
    /// traversed node is a child of i in the position that is "opposite" to pos
    /// (opposite in the sense of opposite_nghbr_position(nghbrPos), see that
    /// function call for details). The traversed positions are saved.
    auto commonParentFound = false;
    auto currentNode = nIdx;
    DBG("starting up-traversal at node: ", nIdx);
    while(!commonParentFound) {
      DBGV((currentNode));
      const auto pIdx = parent(currentNode);
      const auto positionInParent = position_in_parent(currentNode);
      DBGV((parent(currentNode)) (positionInParent)
           (rel_sibling_position(positionInParent,nghbrPos)));

      /// If the parent has a sibling in direction "nghbrPos" we have found
      /// a common parent and are done
      if(is_valid(rel_sibling_position(positionInParent,nghbrPos))) {
        DBG("common parent found -> up-traversal finished");
        commonParentFound = true;
        break;
      }

      if(is_root(pIdx)) {
        /// If the parent is the root and there is no sibling in direction
        /// "nbghrPos", cell has no neighbor and we are done
        DBG("parent is root, and no sibling in direction found!",
            "-> nghbr not found for nIdx: ", nIdx, ", nghbrPos: ", nghbrPos);
        TRACE_OUT();
        return invalid<NodeIdx>();
      } else {
        /// Keep on going up the tree
        traversedPositions.push_back(positionInParent);
        currentNode = pIdx;
      }
    }
    ASSERT(commonParentFound,"No common parent has been found!");

    /// This is the opposite node at the common parent:
    currentNode
        = child(parent(currentNode),
                rel_sibling_position(position_in_parent(currentNode),nghbrPos));

    DBG("starting down traversal at node: ", currentNode);
    /// Traverse the tree back from the opposite node in opposite order to find
    /// the neighbor:
    for(const auto traversedPosition : traversedPositions | boost::adaptors::reversed) {
      const auto nextChildPosition
          = rel_sibling_position(traversedPosition,
                                 opposite_nghbr_position(nghbrPos));
      const auto nextIdx = child(currentNode,nextChildPosition);
      DBGV((currentNode)(traversedPosition)(nextChildPosition)(nextIdx));
      if(is_valid(nextIdx)) {
        currentNode = nextIdx;
      } else {
        DBG("nextIdx: ", nextIdx, " does not exist!",
            " -> nghbr not found for nIdx ", nIdx, ", nghbrPos: " , nghbrPos );
        TRACE_OUT();
        return invalid<NodeIdx>();
      }
    }

    DBG("nghbr found for nIdx: ", nIdx, ", nghbrPos: ", nghbrPos, "!",
        " -> nghbrIdx: " , currentNode);

    assert_valid(currentNode); // Current node must be a valid neighbor
    ASSERT(level(nIdx) == level(currentNode), "Neighbors are not at the same level!");
    TRACE_OUT();
    return currentNode;
  }

  /// \brief Returns the same level neighbors of the node nIdx
  ///
  /// \param [in] nIdx node id of the cell whose neighbors are to be computed.
  /// \returns [IndA<no_samelvl_nghbr_pos()>] Array containing the same level
  /// neighbors of node \p nIdx.
  NodeIdxA<no_samelvl_nghbr_pos()> find_samelvl_nghbrs(const NodeIdx nIdx) const {
    TRACE_IN((nIdx));
    NodeIdxA<no_samelvl_nghbr_pos()> nodeNghbrs;
    for(const auto pos : nghbr_pos()) {
      nodeNghbrs[pos] = find_samelvl_nghbr(nIdx,pos);
    }
    TRACE_OUT();
    return nodeNghbrs;
  }

  ///@}
  //////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////////
  /// \name Node range filters
  ///@{

  /// \brief Returns [RangeFilter] of existing node ids
  inline auto valid() const -> RangeFilter<NodeIdx> {
    return {[&](const NodeIdx nIdx){ return is_valid(nIdx); }};
  }
  /// \brief Returns [RangeFilter] of active node ids
  inline auto active() const -> RangeFilter<NodeIdx> {
    return {[&](const NodeIdx nIdx){ return active_(nIdx); }};
  }
  /// \brief Returns [RangeFilter] of leaf node ids:
  inline auto leafs() const -> RangeFilter<NodeIdx> {
    return {[&](const NodeIdx nIdx){return is_leaf(nIdx);}};
  }
  /// \brief Returns [RangeFilter] of node ids at level "l"
  inline auto atLevel(const SInd l) const -> RangeFilter<NodeIdx> {
    return {[&,l](const NodeIdx nIdx){ return level(nIdx) == l; }};
  }

  inline auto pos_to_nghbr(const NodeIdx nIdx) const -> RangeTransformer<SInd,NodeIdx> {
    return {[&,nIdx](const SInd pos){ return find_samelvl_nghbr(nIdx,pos); }};
  }

  ///@}
  //////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////////
  /// \name Node ranges
  ///@{

  /// \brief Returns [IndRange] of _all_ node position Idxs. That is,
  /// including inactive nodes.
  inline auto all_nodes() const -> Range<NodeIdx> {
    TRACE_IN_();
    TRACE_OUT();
    return {node_begin(), node_end()};
  }
  /// \brief Returns [FilteredRange] of all active node Idxs.
  inline auto nodes() const -> FRange<NodeIdx> {
    TRACE_IN_();
    TRACE_OUT();
    return all_nodes() | active();
  }
  /// \brief Returns [FilteredRange] of all active leaf node Idxs.
  inline auto leaf_nodes() const -> F2Range<NodeIdx> {
    TRACE_IN_();
    TRACE_OUT();
    return nodes() | leafs();
  }
  /// \brief Returns [IndRange] of _all_ child Idxs of node \p nIdx.
  /// That is, including inactive childs.
  inline auto all_childs(const NodeIdx nIdx) const -> Range<NodeIdx*> {
    TRACE_IN((nIdx));
    assert_valid(nIdx);
    TRACE_OUT();
    return boost::make_iterator_range
        (&childrenIds_(nIdx(), 0), &childrenIds_(nIdx(), 0) + no_child_pos());
  }
  /// \brief Returns [FilteredRange] of existing child Idxs of node \p nIdx.
  inline auto childs(const NodeIdx nIdx) const -> FRange<NodeIdx*>{
    TRACE_IN((nIdx));
    assert_valid(nIdx); assert_active(nIdx);
    TRACE_OUT();
    return all_childs(nIdx) | valid();
  }
  /// \brief Returns [IndRange] of all child positions.
  inline auto child_pos() const -> Range<SInd> {
    TRACE_IN_();
    TRACE_OUT();
    return {SInd(0),no_child_pos()}; }

  /// \brief Returns [IndRange] of all nghbr positions.
  inline auto nghbr_pos() const -> Range<SInd> {
    TRACE_IN_();
    TRACE_OUT();
    return {SInd{0},no_samelvl_nghbr_pos()};
  }
  /// \brief Returns [Range] of _all_ nghbr Idxs of node nIdx. That is,
  /// including inactive neighbors.
  inline auto all_samelvl_nghbrs(const NodeIdx nIdx, strict) const
  -> decltype(find_samelvl_nghbrs(nIdx)) {
    TRACE_IN((nIdx));
    assert_valid(nIdx);
    TRACE_OUT();
    return find_samelvl_nghbrs(nIdx);
  }
  inline auto all_samelvl_nghbrs(const NodeIdx nIdx, lazy) const
  -> decltype(nghbr_pos() | pos_to_nghbr(nIdx)){
    TRACE_IN((nIdx));
    assert_valid(nIdx);
    TRACE_OUT();
    return nghbr_pos() | pos_to_nghbr(nIdx);
  }
  template<class T = lazy>
  inline auto all_samelvl_nghbrs(const NodeIdx nIdx, T t = T()) const
  -> decltype(all_samelvl_nghbrs(nIdx,t)) {
    return all_samelvl_nghbrs(nIdx,t);
  }
  /// \brief Returns [FilteredRange] of existing nghbr Idxs of node \p nIdx.
  inline auto samelvl_nghbrs(const NodeIdx nIdx) const
  -> decltype(all_samelvl_nghbrs(nIdx) | valid()) {
    TRACE_IN((nIdx));
    assert_valid(nIdx); assert_active(nIdx);
    TRACE_OUT();
    return all_samelvl_nghbrs(nIdx) | valid();
  }

  ///@}
  //////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////////
  /// \name Solver grids cell ranges
  ///@{

  /// \brief Filters nodes that belong to the \p solverIdx 's grid
  inline auto in_solver(const SolverIdx solverIdx) const -> RangeFilter<NodeIdx> {
    return {[&,solverIdx](const NodeIdx nIdx){ return has_solver(nIdx,solverIdx); }};
  }

  /// \brief Transforms a range of solver ids into a range of cells ids within
  /// each solver's grid for the node \p nIdx
  inline auto solver_pos_to_cell(const NodeIdx nIdx) const
  -> RangeTransformer<SolverIdx,CellIdx> {
    return {[&,nIdx](const SolverIdx solverIdx){ return cell_id(nIdx,solverIdx); }};
  }

  /// \brief Transforms a range of nodes into a range of cell ids belonging to
  /// the \p solverIdx 's grid
  inline auto node_to_cell(const SolverIdx solverIdx) const
  -> RangeTransformer<NodeIdx,CellIdx> {
    return {[&,solverIdx](const NodeIdx nIdx){ return cell_id(nIdx,solverIdx); }};
  }

  /// \brief Returns range of _all_ solver ids.
  inline auto solver_ids() const -> Range<SolverIdx> {
    TRACE_IN_();
    TRACE_OUT();
    return {SolverIdx{0},SolverIdx{solver_capacity()}};
  }

  /// \brief Returns range of _all_ cell ids of node \p nIdx,
  /// including invalid ones!
  inline auto all_cell_ids(const NodeIdx nIdx) const
  -> decltype(solver_ids() | solver_pos_to_cell(nIdx)) {
    TRACE_IN_();
    TRACE_OUT();
    return solver_ids() | solver_pos_to_cell(nIdx);
  }

  /// \brief Return range of all nodes belonging to the \p solverIdx 's grid
  inline auto nodes(const SolverIdx solverIdx) const
  -> decltype(nodes() | in_solver(solverIdx)) {
    return nodes() | in_solver(solverIdx);
  }

  /// \brief Return range of all cells belonging to the \p solverIdx 's grid
  inline auto cells(const SolverIdx solverIdx) const
  -> decltype(nodes(solverIdx) | node_to_cell(solverIdx)) {
    return nodes(solverIdx) | node_to_cell(solverIdx);
  }

  ///@}
  //////////////////////////////////////////////////////////////////////////////


 private:

  //////////////////////////////////////////////////////////////////////////////
  /// \name Container Status Implementation details
  ///@{

  Ind noActiveNodes_;      ///< \f$\#\f$ of active nodes
  Ind maxNoNodes_;         ///< max \f$\#\f$ of nodes
  Ind noNodes_;            ///< \f$\#\f$ of nodes
  NodeIdx lowerFreeNodeBound_; ///< Smallest free node id

  //////////////////////////////////////////////////////////////////////////////
  /// \name Node Data
  ///@{

  /// Alias for data containers:
  template<template <SInd> class V_, SInd nd_ = 1>
  using M = container::Matrix<Implementation,container::matrix::tag::Cell,V_,Ind,SInd,nd_>;

  /// Parent node ids for each node
  M<NodeIdxM>                    parentIds_;

  /// Child node ids for each node
  M<NodeIdxRM,no_child_pos()>    childrenIds_ ;

  /// Indicates if a node is free or in use
  M<BoolMatrix>                  isFree_ ;

  ///@}
  //////////////////////////////////////////////////////////////////////////////

  /// Is the container ready for use, i.e. no modifications are currently
  /// in process and all invariants hold.
  bool isReady_;

  /// \brief Is the \f$\#\f$ of nodes equalto the \f$\#\f$ of active nodes?
  bool is_compact() { return noActiveNodes_ == noNodes_; };

  /// Mapping from node indices to solver cell indices (one for each solver grid)
  EigenDynRowMajor<CellIdx> node2cells_;

  /// \brief Smallest node id that is not in use. If all
  /// node ids are in use, returns size (i.e. one past the end).
  inline NodeIdx free_spot_() const {
    TRACE_IN_();
    const auto spot
        = boost::find_if(Range<NodeIdx>(lowerFreeNodeBound_,node_end()),
                         [&](const NodeIdx i){ return isFree_(i()); })
        - all_nodes().begin();
    TRACE_OUT();
    return NodeIdx{spot};
  }

  /// \brief Writable reference to noActiveNodes_
  inline Ind& size_() { return noActiveNodes_; }
  /// \brief Writable reference no noNodes_
  inline Ind& no_nodes_() { return noNodes_; }

  ///@}
  //////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////////
  /// \name Write access to container variables
  ///@{

  /// \brief Writable reference to parent id of node id \p nIdx
  inline NodeIdx& parent_(const NodeIdx nIdx) { assert_valid(nIdx); return parentIds_(nIdx()); }
  /// \brief Writable reference to child id of \p nIdx at position \p pos
  inline NodeIdx& child_(const NodeIdx nIdx, const SInd pos) {
    TRACE_IN((nIdx)(pos));
    assert_valid(nIdx); assert_child_position(pos);
    TRACE_OUT();
    return childrenIds_(nIdx(),pos);
  }
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
    isFree_(node_begin()()) = false; // activate node before reseting it
    reset_node_(node_begin());
    isFree_(node_begin()()) = false;
    TRACE_OUT();
  }

  /// \brief Resets data of all nodes
  void reset_nodes_(const NodeIdx first, const NodeIdx last) {
    TRACE_IN_();
    for(const auto nIdx : Range<NodeIdx>(first,last)) { reset_node_(nIdx); }
    TRACE_OUT();
  }

  /// \brief Resets all the ids in \p nIdx to invalid<Ind>() and marks node as free.
  void reset_node_(const NodeIdx nIdx) {
    TRACE_IN((nIdx));
    parent_(nIdx) = invalid<NodeIdx>();
    for(const auto pos : child_pos())  { child_(nIdx,pos) = invalid<NodeIdx>(); }
    for(const auto pos : solver_ids()) { cell_id(nIdx,pos) = invalid<CellIdx>(); }
    isFree_(nIdx()) = true;
    TRACE_OUT();
  }

  ///@}
  //////////////////////////////////////////////////////////////////////////////

};


} // hierarchical namespace

template<SInd nd> using Hierarchical = hierarchical::Implementation<nd>;

}} // hom3::container namespace

////////////////////////////////////////////////////////////////////////////////
#undef assert_valid
#undef assert_active
#undef assert_child_position
#undef ENABLE_DBG_
////////////////////////////////////////////////////////////////////////////////
#endif
