#ifndef HOM3_CONTAINER_HIERARCHICAL_IMPLEMENTATION_HPP_
#define HOM3_CONTAINER_HIERARCHICAL_IMPLEMENTATION_HPP_
////////////////////////////////////////////////////////////////////////////////
/// \file \brief Defines the hierarchical container class.
////////////////////////////////////////////////////////////////////////////////
/// Options:
#define ENABLE_DBG_ 0
#include "../../misc/dbg.hpp"
////////////////////////////////////////////////////////////////////////////////

namespace container { namespace hierarchical {

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

/// Stencil of relative sibling positions
static constexpr std::array<std::array<SInd,6>,8>  rel_sibling_position_arr {{
    //    0       1       2        3         4        5   | nghbrPos
    {{ iSInd(),       1, iSInd(),       2, iSInd(),       4 }}, // child 0
    {{       0, iSInd(), iSInd(),       3, iSInd(),       5 }}, // child 1
    {{ iSInd(),       3,       0, iSInd(), iSInd(),       6 }}, // child 2
    {{       2, iSInd(),       1, iSInd(), iSInd(),       7 }}, // child 3
    {{ iSInd(),       5, iSInd(),       6,       0, iSInd() }}, // child 4
    {{       4, iSInd(), iSInd(),       7,       1, iSInd() }}, // child 5
    {{ iSInd(),       7,       4, iSInd(),       2, iSInd() }}, // child 6
    {{       6, iSInd(),       5, iSInd(),       3, iSInd() }}  // child 7
  }};

/// Stencil of opposite neighbor positions
static constexpr std::array<SInd,6> opposite_nghbr_position_arr{{ 1, 0, 3, 2, 5, 4 }};


/// \brief Hierarchical Node Container: "nd"-tree that handles
/// "nc" container ids per node.
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
/// Conventions:
///
/// Identifiers:
/// pId -> parentId
/// cId -> childId / cellId (in container)
/// cPos -> child position (in parent)
/// nId -> nodeId
/// nPos -> neighbor position (relative to node)
/// eId -> endId
///
/// p -> Predicate/Kernel
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

  /// \brief Constructs a grid container that can maximally hold "nn" nodes.
  ///
  /// \param [in] nn maximum number of nodes that the container can hold.
  Implementation(const Ind nn, const SInd nc)
      : noActiveNodes_{0}, maxNoNodes_{nn}, lowerFreeNodeBound_{0},
    parents_(this,"parents"), firstChild_(this,"firstChild"),
    hasChildInPos_(this,"hasChildInPos"), isFree_(this,"isFree"),
    isReady_{false}, g2l_{nn,nc} {
    TRACE_IN((nn)(nc));

    initialize_(parents_,childs_,isFree_);

    reset_all_nodes_();

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
  /// \name Global-to-local mapping
  ///@{

  /// \brief Local id of node \p globalId within \p container
  inline Ind local_id(const Ind globalId, const SInd container) const {
    assert_valid(globalId); assert_active(globalId);
    return g2l_(globalId,container);
  }
  /// \brief Reference to the local id of node \p globalId within \p container
  inline Ind& local_id(const Ind globalId, const SInd container) {
    assert_valid(globalId); assert_active(globalId);
    return g2l_(globalId,container);
  }
  /// \brief Does \p container have a localId of node \p globalId ?
  inline bool has_container(const Ind globalId, const SInd container) const {
    assert_valid(globalId); assert_active(globalId);
    return is_valid(g2l_(globalId,container));
  }
  /// \brief \f$\#\f$ of containers that can be registered in the grid
  inline SInd container_capacity() const {
    return g2l_.cols();
  }

  ///@}
  //////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////////
  /// \name Node operations
  ///@{

  /// \brief Is node "nId" an active node?
  /// (active means "in use", i.e. not free)
  inline bool active_(const Ind nId) const {
    TRACE_IN((nId));
    assert_valid(nId);
    TRACE_OUT();
    return !isFree_(nId);
  }
  /// \brief Is the node a root node?
  /// (right now there is only a single root node)
  inline bool is_root(const Ind nId) const {
    TRACE_IN((nId));
    assert_valid(nId);
    TRACE_OUT();
    return !is_valid(parent(nId));
  }
  /// \brief Is node "nId" a leaf node?
  inline bool is_leaf(const Ind nId) const {
    TRACE_IN((nId));
    assert_valid(nId); assert_active(nId);
    TRACE_OUT();
    return no_childs(nId) == 0;
  }
  /// \brief \f$\#\f$ of leaf nodes in the container
  inline Ind no_leaf_nodes() const { TRACE_IN_(); TRACE_OUT();
    return boost::distance(leaf_nodes());
  }

  /// \brief Id of parent of "nId"
  inline const Ind& parent(const Ind nId) const { assert_valid(nId); return parents_(nId); }

  /// \brief Id of the child of "nId" locatet at position "pos"
  inline const Ind& child(const Ind nId, const SInd pos) const {
    TRACE_IN((nId)(pos));
    assert_valid(nId); assert_child_position(pos);
    TRACE_OUT();
    if(!has_child_at_pos(nId,pos)) { return iInd(); }
    return firstChild_(nId)
        + boost::count_if(boost::counting_range(0,pos),[&](const SInd p){
            return has_child_at_pos(nId,p);
          });
  }

  /// \brief \f$\#\f$ of children of "nId"
  inline SInd no_childs(const Ind nId) const {
    TRACE_IN((nId));
    assert_valid(nId);
    TRACE_OUT();
    return boost::count(all_child_pos(),[&](const SInd pos){
        return has_child_at_pos(nId,pos);
      });
  }
  inline bool has_child_at_pos(const Ind nId, const SInd pos) const {
    return hasChildAtPos_(nId,pos);
  }
  /// \brief Does "nId" have children?
  inline bool has_children(const Ind nId) const { return !is_leaf(nId); }
  /// \brief Does "nId" have _all_ children? I.e.
  /// is \f$\#\f$ of children of "nId" == \f$\#\f$ of children positions?
  inline bool has_all_children(const Ind nId) const {
    TRACE_IN((nId));
    assert_valid(nId);
    TRACE_OUT();
    return no_childs(nId) == no_child_pos();
  }
  /// \brief Position in parent of node "nId"
  SInd position_in_parent(const Ind nId) const {
    TRACE_IN((nId));
    assert_valid(nId); assert_active(nId);
    const auto pId = parent(nId);
    assert_valid(pId); assert_active(pId);
    const SInd res = boost::find_if(child_pos(), [&](const SInd cPos){
        return child(pId,cPos) == nId; }) - child_pos().begin();
    TRACE_OUT();
    return res;
  }

  /// \brief #of same level neighbors of node \p nId
  inline SInd no_samelvl_nghbrs(const Ind nId) const {
    TRACE_IN((nId));
    assert_valid(nId);
    TRACE_OUT();
    return boost::distance(samelvl_nghbrs(nId));
  }

  /// \brief Opposite nghbr position. I.e. the neighbor with
  /// position "pos" of a given nId, has nId as neighbor in its opposite
  /// position.
  ///
  /// \warning returns iSInd if out of bounds. Bound check
  /// not performed here.
  static constexpr SInd opposite_nghbr_position(const SInd pos) {
    return DBG_EXPR(pos < 6 ?)
        opposite_nghbr_position_arr[pos]
        DBG_EXPR(: iSInd());
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
           DBG_EXPR(: iSInd());
  }
  ///@}
  //////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////////
  /// \name Modifying algorithms
  ///@{

  /// \brief Inserts a node into parent "pId" at position "pos".
  ///
  /// \complexity O(1) if the tree is_compact()
  /// \complexity average O(1) if the tree is not Compact
  Ind insert_node(const Ind pId, const SInd pos) {
    TRACE_IN((pId)(pos));
    assert_active(pId);
    assert_valid(!child(pId,pos)); // position must be free!
    ASSERT(level(pId) + 1 < max_no_levels(), "can't refine > max_no_levels()!");
    const auto nId = is_compact() ? size() : free_spot_();
    ++no_nodes_();
    ++lowerFreeNodeBound_;
    if(nId == size()) {
      ++size_();
    }
    reset_node_(nId);
    child_(pId,pos) = nId;
    parent_(nId) = pId;
    isFree_(nId) = false;
    TRACE_OUT();
    return nId;
  }

  /// \brief Removes node "nId" from its parent.
  ///
  /// \complexity: O(1)
  Ind remove_node(const Ind nId) {
    TRACE_IN((nId));
    assert_active(nId);
    const auto pId = parent(nId);
    ASSERT(is_valid(pId), "Node has no parent!");
    const auto pos = position_in_parent(nId);
    ASSERT(is_valid(pos), "Node is not child of parent!");
    child_(pId,pos) = iInd();
    reset_node_(nId);
    --no_nodes_();
    if(nId == size()) {
      --size();
    }
    isReady_ = false;
    lowerFreeNodeBound_ = std::min(lowerFreeNodeBound_,nId);
    TRACE_OUT();
    return pId;
  }

  /// \brief Refines a node isotropically, i.e. into 4 (in 2D) or 8 (in 3D)
  /// children. It optionally takes a predicate "p(posInParent)" that is
  /// executed after each child is inserted.
  ///@{
  template<class P>
  void refine_node(const Ind nId, P&& p) {
    TRACE_IN((nId));
    for(auto pos : child_pos()) {
      insert_node(nId,pos); // note insert_node returns nId of child
      p(pos);               // but computing it from the nId is easy with pos.
    }
    isReady_ = false;
    TRACE_OUT();
  }
  void refine_node(const Ind nId) { refine_node(nId,[](const SInd){}); }
  ///@}

  ///@}
  //////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////////
  /// \name Non-modifying algorithms
  ///@{

  /// \brief Computes the level of the node "nId"
  /// \complexity
  ///   - O(\f$\#\f$l): linear in the number of levels,
  ///   - O(log(\f$\#\f$size)): logarithmic in the number of cells
  inline SInd level(Ind nId) const {
    TRACE_IN((nId));
    assert_active(nId);
    SInd level_ = 0;
    nId = parent(nId);
    while(is_valid(nId)) {
      assert_active(nId);
      nId = parent(nId);
      ++level_;
    }
    ASSERT(level_ < max_no_levels(), "error!");
    TRACE_OUT();
    return level_;
  }

  /// \brief Computes the same level neighbor of the node "nId" located at
  /// its relative position "pos".
  ///
  /// \param [in] nId  node id of the cell whose neighbor is to be determined.
  /// \param [in] pos  position of the required neighbor w.r.t. the node nId.
  ///
  /// \returns [Ind] node id of the neighbor located at pos w.r.t. the node nId.
  /// If no neighbor is found it returns iInd().
  ///
  /// \complexity unknown (worst case is probably O(\f$\#\f$l) = 2 * \f$\#\f$l,
  /// need to traverse the tree up to the root and back down again).
  /// \todo Improve complexity and profile!
  ///
  /// \requires a balanced tree that satisfies 2:1 rule.
  Ind find_samelvl_nghbr(const Ind nId, const SInd nghbrPos) const {
    TRACE_IN((nId)(nghbrPos));
    DBG("start samelvl_nghbr | nId: ", nId, " | nghbrPos ", nghbrPos);

    /// The root cell has no neighbors:
    if(is_root(nId)) {
      DBG("nId ", nId, " is a root cell | no nghbr found in nghbrPos: ", nghbrPos);
      TRACE_OUT();
      return iInd();
    }

    // Create vector of reserved positions on the stack
    memory::stack::arena<SInd,max_no_levels()> stackMemory;
    auto traversedPositions = memory::stack::make<std::vector>(stackMemory);
    traversedPositions.reserve(max_no_levels() /*level(nId)*/);

    /// Travel up the tree until a node i is found, s.t. the previously
    /// traversed node is a child of i in the position that is "opposite" to pos
    /// (opposite in the sense of opposite_nghbr_position(nghbrPos), see that
    /// function call for details). The traversed positions are saved.
    auto commonParentFound = false;
    auto currentNode = nId;
    DBG("starting up-traversal at node: ", nId);
    while(!commonParentFound) {
      DBGV((currentNode));
      const auto pId = parent(currentNode);
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

      if(is_root(pId)) {
        /// If the parent is the root and there is no sibling in direction
        /// "nbghrPos", cell has no neighbor and we are done
        DBG("parent is root, and no sibling in direction found!",
            "-> nghbr not found for nId: ", nId, ", nghbrPos: ", nghbrPos);
        TRACE_OUT();
        return iInd();
      } else {
        /// Keep on going up the tree
        traversedPositions.push_back(positionInParent);
        currentNode = pId;
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
      const auto nextId = child(currentNode,nextChildPosition);
      DBGV((currentNode)(traversedPosition)(nextChildPosition)(nextId));
      if(is_valid(nextId)) {
        currentNode = nextId;
      } else {
        DBG("nextId: ", nextId, " does not exist!",
            " -> nghbr not found for nId ", nId, ", nghbrPos: " , nghbrPos );
        TRACE_OUT();
        return iInd();
      }
    }

    DBG("nghbr found for nId: ", nId, ", nghbrPos: ", nghbrPos, "!",
        " -> nghbrId: " , currentNode);

    assert_valid(currentNode); // Current node must be a valid neighbor
    ASSERT(level(nId) == level(currentNode), "Neighbors are not at the same level!");
    TRACE_OUT();
    return currentNode;
  }

  /// \brief Returns the same level neighbors of the node nId
  ///
  /// \param [in] nId node id of the cell whose neighbors are to be computed.
  /// \returns [IndA<no_samelvl_nghbr_pos()>] Array containing the same level
  /// neighbors of node "nId".
  IndA<no_samelvl_nghbr_pos()> find_samelvl_nghbrs(const Ind nId) const {
    TRACE_IN((nId));
    IndA<no_samelvl_nghbr_pos()> nodeNghbrs;
    for(const auto pos : nghbr_pos()) {
      nodeNghbrs[pos] = find_samelvl_nghbr(nId,pos);
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
  inline auto exist() const -> RangeFilter<Ind> {
    return {[&](const Ind nId){ return is_valid(nId); }};
  }
  /// \brief Returns [RangeFilter] of active node ids
  inline auto active() const -> RangeFilter<Ind> {
    return {[&](const Ind nId){ return active_(nId); }};
  }
  /// \brief Returns [RangeFilter] of leaf node ids:
  inline auto leafs() const -> RangeFilter<Ind> {
    return {[&](const Ind nId){return is_leaf(nId);}};
  }
  /// \brief Returns [RangeFilter] of node ids at level "l"
  inline auto atLevel(const SInd l) const -> RangeFilter<Ind> {
    return {[&](const Ind nId){ return level(nId) == l; }};
  }

  inline auto pos_to_nghbr(const Ind nId) const -> RangeTransformer<SInd,Ind> {
    return {[&](const SInd pos){ return find_samelvl_nghbr(nId,pos); }};
  }

  ///@}
  //////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////////
  /// \name Node ranges
  ///@{

  /// \brief Returns [IndRange] of _all_ node position Ids. That is,
  /// including inactive nodes.
  inline auto all_nodes() const -> Range<Ind> {
    TRACE_IN_();
    TRACE_OUT();
    return {Ind(0),size()};
  }
  /// \brief Returns [FilteredRange] of all active node Ids.
  inline auto nodes() const -> FRange<Ind> {
    TRACE_IN_();
    TRACE_OUT();
    return all_nodes() | active();
  }
  /// \brief Returns [FilteredRange] of all active leaf node Ids.
  inline auto leaf_nodes() const -> F2Range<Ind> {
    TRACE_IN_();
    TRACE_OUT();
    return nodes() | leafs();
  }
  /// \brief Returns [IndRange] of _all_ child Ids of node "nId".
  /// That is, including inactive childs.
  inline auto all_childs(const Ind nId) const -> Range<Ind*> {
    TRACE_IN((nId));
    assert_valid(nId);
    TRACE_OUT();
    return boost::make_iterator_range
        (&childs_(nId, 0), &childs_(nId, 0) + no_child_pos());
  }
  /// \brief Returns [FilteredRange] of existing child Ids of node "nId".
  inline auto childs(const Ind nId) const -> FRange<Ind*>{
    TRACE_IN((nId));
    assert_valid(nId); assert_active(nId);
    TRACE_OUT();
    return all_childs(nId) | exist();
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
  /// \brief Returns [Range] of _all_ nghbr Ids of node nId. That is,
  /// including inactive neighbors.
  inline auto all_samelvl_nghbrs(const Ind nId, strict) const
  -> decltype(find_samelvl_nghbrs(nId)) {
    TRACE_IN((nId));
    assert_valid(nId);
    TRACE_OUT();
    return find_samelvl_nghbrs(nId);
  }
  inline auto all_samelvl_nghbrs(const Ind nId, lazy) const
  -> decltype(nghbr_pos() | pos_to_nghbr(nId)){
    TRACE_IN((nId));
    assert_valid(nId);
    TRACE_OUT();
    return nghbr_pos() | pos_to_nghbr(nId);
  }
  template<class T = lazy>
  inline auto all_samelvl_nghbrs(const Ind nId, T t = T()) const
  -> decltype(all_samelvl_nghbrs(nId,t)) {
    return all_samelvl_nghbrs(nId,t);
  }
  /// \brief Returns [FilteredRange] of existing nghbr Ids of node "nId".
  inline auto samelvl_nghbrs(const Ind nId) const -> decltype(all_samelvl_nghbrs(nId) | exist()) {
    TRACE_IN((nId));
    assert_valid(nId); assert_active(nId);
    TRACE_OUT();
    return all_samelvl_nghbrs(nId) | exist();
  }

  ///@}
  //////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////////
  /// \name Container ranges
  ///@{

  inline auto is_valid_container() const -> RangeFilter<SInd> {
    return {[](const SInd container){ return is_valid(container); }};
  }

  inline auto global_to_local(const Ind cId) const -> RangeTransformer<SInd,Ind> {
    return {[&](const SInd pos){ return g2l_(cId,pos); }};
  }

  /// \brief Returns range of _all_ container positions.
  inline auto container_pos() const -> Range<SInd> {
    TRACE_IN_();
    TRACE_OUT();
    return {SInd(0),container_capacity()};
  }

  /// \brief Returns range of all containers of cId, including inactive ones!
  // inline auto all_containers(const Ind cId) const -> FRange<SInd> {
  inline auto all_containers(const Ind cId) const
      -> decltype(container_pos() | global_to_local(cId)) {
    TRACE_IN_();
    TRACE_OUT();
    return container_pos() | global_to_local(cId);
  }

  /// \brief Returns range of all active containers of cIdnode Ids.
  inline auto containers(const Ind cId) const
      -> decltype(all_containers(cId) | is_valid_container()) {
    TRACE_IN_();
    TRACE_OUT();
    return all_containers(cId) | is_valid_container();
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
  Ind lowerFreeNodeBound_; ///< Smallest free node id

  //////////////////////////////////////////////////////////////////////////////
  /// \name Node Data
  ///@{

  /// Alias for data containers:
  template<template <SInd> class V_, SInd nd_ = 1>
  using M = container::Matrix<Implementation,container::matrix::tag::Cell,V_,nd_>;

  /// Parent node ids for each node
  M<IndM> parents_;

  /// First Child node id for each node
  M<IndM> firstChild_;

  /// Has child in pos?
  M<BoolMatrix,2*nd> hasChildInPos_;
  
  /// Indicates if a node is free or in use
  M<BoolMatrix> isFree_ ;

  ///@}
  //////////////////////////////////////////////////////////////////////////////

  /// Is the container ready for use, i.e. no modifications are currently
  /// in process and all invariants hold.
  bool isReady_;

  /// \brief Is the \f$\#\f$ of nodes equalto the \f$\#\f$ of active nodes?
  bool is_compact() { return noActiveNodes_ == noNodes_; };

  /// Global to local mapping (for each container)
  EigenDynRowMajor<Ind> g2l_;

  /// \brief Smallest node id that is not in use. If all
  /// node ids are in use, returns size (i.e. one past the end).
  inline Ind free_spot_() const {
    TRACE_IN_();
    const auto spot
        = boost::find_if(Range<Ind>(lowerFreeNodeBound_,size()),
                         [&](const Ind i){ return isFree_(i); })
        - all_nodes().begin();
    TRACE_OUT();
    return spot;
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

  /// \brief Writable reference to parent id of node id "nId"
  inline Ind& parent_(const Ind nId) { assert_valid(nId); return parents_(nId); }
  /// \brief Writable reference to child id of "nId" at position "pos"
  inline Ind& child_(const Ind nId, const SInd pos) {
    TRACE_IN((nId)(pos));
    assert_valid(nId); assert_child_position(pos);
    TRACE_OUT();
    return childs_(nId,pos);
  }
  ///@}
  //////////////////////////////////////////////////////////////////////////////


  //////////////////////////////////////////////////////////////////////////////
  /// \name Initialization routines
  ///@{

  /// \brief Calls init method for containers of variables
  ///
  /// The init method sets ownership by setting the container pointer to this
  template<class H> void initialize_(H&& head) { head.init(this); }
  template<class H, class... Ts> void initialize_(H&& head, Ts&&... tail) {
    initialize_(std::forward<H>(head));
    initialize_(std::forward<Ts>(tail)...);
  }

  /// \brief Creates an initial root node
  ///
  /// \warning The container must be empty!
  void initialize_root_node_() {
    TRACE_IN_();
    ASSERT(empty(), "Container is not empty!");
    ++size_(); ++no_nodes_();
    reset_node_(0);
    isFree_(0) = false;
    TRACE_OUT();
  }


  /// \brief Resets data of all nodes
  void reset_all_nodes_() {
    TRACE_IN_();
    for(const auto nId : all_nodes()) { reset_node_(nId); }
    TRACE_OUT();
  }

  /// \brief Resets all the ids in "nId" to iInd() and marks node as free.
  void reset_node_(const Ind nId) {
    TRACE_IN((nId));
    parent_(nId) = iInd();
    isFree_(nId) = true;
    for(const auto pos : child_pos()) { child_(nId,pos) = iInd(); }
    for(const auto pos : container_pos()) { g2l_(nId,pos) = iInd(); }
    TRACE_OUT();
  }

  ///@}
  //////////////////////////////////////////////////////////////////////////////

};


} // hierarchical namespace

template<SInd nd> using Hierarchical = hierarchical::Implementation<nd>;

} // container namespace

#undef ENABLE_DBG_
////////////////////////////////////////////////////////////////////////////////
#endif
