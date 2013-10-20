#ifndef HOM3_CONTAINER_SEQUENTIAL_IMPLEMENTATION_HPP_
#define HOM3_CONTAINER_SEQUENTIAL_IMPLEMENTATION_HPP_
////////////////////////////////////////////////////////////////////////////////
/// \file \brief Implementation of the sequential container
////////////////////////////////////////////////////////////////////////////////
/// Includes:
#include <type_traits>
#include <limits>
#include <algorithm>
////////////////////////////////////////////////////////////////////////////////
#include "containers/sequential/traits.hpp"
#include "containers/sequential/iterator.hpp"
#include "containers/sequential/algorithm.hpp"
////////////////////////////////////////////////////////////////////////////////
/// Options:
#define ENABLE_DBG_ 0
#include "misc/dbg.hpp"
////////////////////////////////////////////////////////////////////////////////
/// File macros:
#define assert_variable_node_container()                                  \
  static_assert(std::is_same<container_type, tag::variable_nodes>::value, \
                "Function defined only for variable_node_size containers!")

#define assert_in_cell_range(cIdx)                                      \
  ASSERT(in_cell_range(cIdx),                                           \
         "Cell index " << cIdx                                          \
         << " is out of bounds [" << first() << ","                     \
                                  << size() <<  ").")

#define assert_in_node_range(nIdx)                                      \
  ASSERT(in_node_range(nIdx),                                           \
         "Node index " << nIdx                                          \
         << " is out of bounds [" << first_node() << ","                \
                                  << last_node() <<  ").")

#define assert_in_cell_node_range(cIdx, nIdx)                           \
  ASSERT(in_cell_node_range(cIdx, nIdx),                                \
         "Node index " << nIdx                                          \
         << " is out of cIdx: " << cIdx                                 \
         << " 's bounds [" << first_node(cIdx) << ","                   \
                           << last_node(cIdx) <<  ").")

#define assert_valid_node_range(fromIdx, toIdx) \
  ASSERT(is_valid_node_range(fromIdx, toIdx),   \
         "Invalid node range [" << fromIdx      \
         << ", " << toIdx << ").")
////////////////////////////////////////////////////////////////////////////////

namespace hom3 { namespace container {

/// \brief Sequential container
namespace sequential {


/// \brief Sequential container (base implementation)
///
/// Defines the interface of cell containers and implements
/// common functionality.
///
/// The pourpose of this class is to reduce the boilerplate
/// required to implement cell containers
///
template<class Cells> struct Implementation {
  /// Traits:
  using iterator       = Iterator<Cells>;
  using value_type     = typename traits<Cells>::value_type;
  using reference      = typename traits<Cells>::reference;
  using const_iterator = const iterator;
  using container_type = typename traits<Cells>::container_type;
  using cell_size_type = Ind;
  using node_size_type = Ind;
  using CIdx           = typename traits<Cells>::cell_index_type;
  using NIdx           = decltype(tti::node_index_type<traits<Cells>, Ind>());


  /// \name Constructors
  ///@{
  Implementation(const cell_size_type n,  const String name) noexcept
    : Implementation(n, 0, name) {}
  Implementation(const cell_size_type ne, const node_size_type nn,
                 const String name) noexcept
    : maxCellSize_(ne)
    , maxNodeSize_(nn)
    , cellSize_(0)
    , nodeSize_(0)
    , nodes_(c(), "nodes")
    , name_(name) {
    TRACE_IN((ne)(nn));

    if (capacity() == 0) {
      TERMINATE("Empty container!");
    }

    if (std::is_same<container_type, tag::variable_nodes>::value
        && node_capacity() < capacity()) {
      TERMINATE("Less nodes than cells!");
    }

    TRACE_OUT();
  }

  Implementation(const Implementation& other) noexcept
    : maxCellSize_(other.capacity())
    , maxNodeSize_(other.node_capacity())
    , cellSize_(other.size())
    , nodeSize_(other.node_size())
    , nodes_(other.nodes_)
    , name_(other.name_) {
    TRACE_IN_();
    init_nodes_(container_type());
    TRACE_OUT();
  }

  Implementation(Implementation&& other) noexcept
    : maxCellSize_(other.capacity())
    , maxNodeSize_(other.node_capacity())
    , cellSize_(other.size())
    , nodeSize_(other.node_size())
    , nodes_(other.nodes_) {
    TRACE_IN_();
    TRACE_OUT();
  }

  /// \brief Initializes a single column type
  template<class Variable> void initialize(Variable&& variable) noexcept
  { variable.init(c()); }
  /// \brief Initialize the column types by calling init for each member  (each
  /// member's init function takes the container pointer and uses it to check
  /// ownership)
  template<class Variable, class... Variables>
  void initialize(Variable&& variable, Variables&&... variables) noexcept {
    initialize(std::forward<Variable>(variable));
    initialize(std::forward<Variables>(variables)...);
  }
  ///@}

  /// \name Size functions
  ///@{

  /// \brief Returns the maximum #of cells that can be stored in the container
  inline cell_size_type capacity() const noexcept { return maxCellSize_; }
  /// \brief Returns the maximum #of nodes that can be stored in the container
  inline node_size_type node_capacity() const noexcept { return maxNodeSize_; }
  /// \brief Is the container empty (i.e. does it have no cells?)
  inline bool empty() const noexcept { return cellSize_ == 0; }
  /// \brief Returns the current #of cells stored in the container
  inline cell_size_type size() const noexcept { return cellSize_; }
  /// \brief Returns the current #of nodes stored in the container
  inline node_size_type node_size() const noexcept { return nodeSize_; }
  /// \brief Returns the current #of nodes of the element \p cIdx
  inline node_size_type node_size(const CIdx cIdx) const noexcept {
    assert_variable_node_container();
    ASSERT(first_node(cIdx) <= last_node(cIdx), "Invalid cell node range!");
    return last_node(cIdx) - first_node(cIdx);
  }
  ///@}

  /// \name Cell iterators
  ///@{

  /// \brief Returns a cell iterator to the first container element
  inline iterator        begin()       noexcept { return {c(), first()}; }
  /// \brief Returns a const cell iterator to the first container element
  inline const_iterator cbegin() const noexcept { return {c(), first()}; }

  /// \brief Returns a cell iterator to one-past-the-last container element
  inline iterator          end()       noexcept { return {c(), last()}; }
  /// \brief Returns a const cell iterator to one-past-the-last element
  inline const_iterator   cend() const noexcept { return {c(), last()}; }

  ///@}

  /// \name Indexed cell access
  ///@{

  /// \brief Returns index of the first container element (always 0)
  inline CIdx front() const noexcept { return CIdx{0}; }
  /// \brief Returns index of the last container element (size()-1)
  inline CIdx back()  const noexcept { return CIdx{size() - 1}; }
  /// \brief Returns index of the first container element (always 0)
  inline CIdx first() const noexcept { return CIdx{0}; }
  /// \brief Returns index of the one-past-end container element (size())
  inline CIdx last()  const noexcept { return CIdx{size()}; }
  /// \brief Random access: returns a reference to the element \p i
  inline reference operator[](const CIdx i) noexcept {
    assert_in_cell_range(i);
    return {c(), i};
  }
  inline reference operator[](const CIdx i) const noexcept {
    assert_in_cell_range(i);
    return {c(), i};
  }
  /// \brief Is cell index \p i in range?
  inline bool in_cell_range(const CIdx i) const noexcept
  { return primitive_cast(i) < size(); }
  ///@}

  /// \name Indexed ranges
  ///@{

  /// \brief Returns [IndRange] of _all_ container cell Ids.
  inline auto all_cells() const noexcept -> Range<CIdx>
  { return {first(), last()}; }
  ///@}

  ///@}

  /// \name Indexed node access
  ///@{

  /// \brief Returns index of the first node of the cell \p cIdx
  inline NIdx first_node(const CIdx cIdx = 0) const noexcept
  { assert_variable_node_container(); return nodes_(cIdx); }

  /// \brief Returns index of the last node of the cell \p cIdx
  ///
  /// Note: In the last cell "last_node(size()-1)" calls "nodes_(size())".
  /// Although the cell with cIdx == size() does not exists, the container has
  /// an extra _magic_ node that marks the end of the last cells' nodes.
  ///
  ///@{
  inline NIdx last_node(const CIdx cIdx) const noexcept
  { assert_variable_node_container(); return nodes_(cIdx + 1);    }
  inline NIdx last_node()                const noexcept
  { assert_variable_node_container(); return last_node(size() - 1); }
  ///@}

  /// \brief Is node \p i in node range?
  inline bool in_node_range(const NIdx i) const noexcept
  { return i < last_node(); }
  /// \brief Is node \p nIdx in the node range of cell \p cIdx?
  inline bool in_cell_node_range
  (const CIdx cIdx, const NIdx nIdx) const noexcept
  { return nIdx >= first_node(cIdx) && nIdx < last_node(cIdx); }
  /// \brief Is node range [\p from, \p toIdx) valid?
  inline bool is_valid_node_range
  (const NIdx fromIdx, const NIdx toIdx) const noexcept
  { return fromIdx < toIdx; }
  ///@}

  /// \name Append/delete cells
  ///@{

  /// \brief Inserts \p noCells with \p noNodesPerCell at the container end
  /// and resets their variables
  inline CIdx push_cell(const cell_size_type noCells = 1,
                        const node_size_type noNodesPerCell = 1) noexcept {
    auto first = push_cell_(noCells, noNodesPerCell, container_type());
    reset_cells_(last() - CIdx{noCells}, last());
    return first;
  }
  /// \brief Removes last \p i cells from the container
  inline void pop_cell(const cell_size_type i = 1) noexcept
  { pop_cell_(i, container_type()); }
  ///@}

  /// \name Append/delete nodes
  ///@{

  /// \brief Inserts \p i nodes at the end of the container
  void push_node(const node_size_type i = 1) noexcept {
    assert_variable_node_container();
    node_size_() += i;
    last_node_() += i;
    ASSERT(node_size() <= node_capacity(), "Container out of memory.");
  }
  /// \brief Removes \p i nodes at the end of the container
  void pop_node(const node_size_type i = 1) noexcept {
    assert_variable_node_container();
    node_size_() -= i;
    last_node_() -= i;
    ASSERT(node_size() >= size(), "Less nodes than cells!");
  }
  ///@}

  /// \name Algorithms
  ///@{

  /// \brief Removes elements in range \p c that satisfy the predicate
  /// \p p
  ///
  /// \algorithm mutating
  /// \complexity O(n) where n = #of elements in the range \p c
  template<class Predicate>
  auto erase_remove_if(Predicate&& p) noexcept -> CIdx {
    const auto last_ = last();
    auto first = algorithm::find_if(CIdx{0}, last_, p);

    if (first == last_) { return first; }

    auto next = first; ++first;
    for (; first != last_; ++first) {
      if (!p(first)) {
        copy_cell(first, next);
        ++next;
      }
    }

    pop_cell(last_ - next);
    return next;
  }

  /// \brief Copies the element at \p fromCIdx to position \p toCIdx
  ///
  /// \complexity O(1) for fixed_nodecontainers
  /// \complexity O(N) for variable_nodes containers
  void copy_cell(const CIdx fromCIdx, const CIdx toCIdx) noexcept {
    ASSERT(fromCIdx != toCIdx, "Suspicious | Untested behaviour");
    copy_cell_(fromCIdx, toCIdx, container_type());
  }

  ///@}

  /// \name Extra functionality
  ///@{

  /// \brief Swap containers
  /// \todo should swap be limited to containers of equal maxCell/Node sizes?
  friend void swap(Cells& first, Cells& second) noexcept {
    using std::swap;
    TRACE_IN_();

    /// Swaps container variables:
    swap(first.maxCellSize_, second.maxCellSize_);
    swap(first.maxNodeSize_, second.maxNodeSize_);
    swap(first.cellSize_   , second.cellSize_);
    swap(first.nodeSize_   , second.nodeSize_);
    swap(first.nodes_      , second.nodes_);

    /// User defined function for swapping user-defined variables:
    swap_containers(first, second);

    TRACE_OUT();
  }
  ///@}

  String name() const noexcept { return name_; }

 private:
  // Note: max sizes aren't const to allow swapping containers of different
  // sizes \todo remove swapping of containers and make them const
  cell_size_type maxCellSize_;  ///< max #of cells + 1
  node_size_type maxNodeSize_;  ///< max #of nodes + 1

  cell_size_type cellSize_;  ///< #of cells + 1
  node_size_type nodeSize_;  ///< #of nodes + 1

  /// Container of node indicies
  NodeIndices<Cells, CIdx> nodes_;

  const String name_;

  /// \brief Returns a pointer to the underlying container (CRTP)
  inline Cells* c() noexcept { return static_cast<Cells*>(this); }

  /// \brief Generic initialization of both fixed and variable node containers
  ///@{
  void init_nodes_(tag::fixed_nodes)    noexcept {                     }
  void init_nodes_(tag::variable_nodes) noexcept { initialize(nodes_); }
  ///@}

  /// \name Allow to modify the container private variables
  ///@{
  inline cell_size_type& size_() noexcept { return cellSize_; }
  inline node_size_type& node_size_() noexcept
  { assert_variable_node_container(); return nodeSize_; }
  inline NIdx& first_node_(const CIdx cIdx = 0) noexcept
  { assert_variable_node_container(); return nodes_(cIdx);     }
  inline NIdx& last_node_(const CIdx cIdx) noexcept
  { assert_variable_node_container(); return nodes_(cIdx + 1); }
  inline NIdx& last_node_() noexcept
  { assert_variable_node_container(); return last_node_(size() - 1); }
  ///@}

  /// Finds the cIdx of a given node (slow: for debugging purposes only!)
  inline CIdx find_cell_with_node(const NIdx nIdx) const noexcept {
    assert_variable_node_container();
    for (auto cId : all_cells()) {
      if (first_node(cId) <= nIdx && last_node(cId) > nIdx) {
        return cId;
      }
    }
    return size();
  }

  /// \name Implementation details of append/delete functions
  ///@{
  inline CIdx push_cell_(const cell_size_type noCells,
                         node_size_type, tag::fixed_nodes) noexcept {
    const auto first = last();
    size_() += noCells;
    ASSERT(size() <= capacity(), "Container out of memory.");
    return first;
  }

  inline CIdx push_cell_
  (const cell_size_type noCells, const node_size_type noNodesPerCell,
  tag::variable_nodes) noexcept {
    const auto first = last();
    auto cIdx = first;
    auto nIdx = node_size();
    size_() += noCells;
    node_size_() += noCells * noNodesPerCell;
    for (; cIdx < size(); ++cIdx) {
      first_node_(cIdx) = nIdx;
      nIdx += noNodesPerCell;
    }
    last_node_() = node_size();
    ASSERT(size() <= capacity(), "Container out of memory.");
    return first;
  }

  void pop_cell_(const CIdx i, tag::fixed_nodes) noexcept {
    ASSERT(!empty(), "Container is already empty!");
    ASSERT(size() - i > 0, "Container doesn't have enough cells to pop!");
    size_() -= i;
  }

  void pop_cell_(const CIdx i, tag::variable_nodes) noexcept {
    ASSERT(!empty(), "Container is already empty!");
    ASSERT(size() - i > 0, "Container doesn't have enough cells to pop!");
    last_node_() = last_node(size() - i - 1);
    node_size_() = last_node();
    size_() -= i;
  }

  /// \brief Shifts node range ["fromNIdx","toNIdx") up "steps" times
  ///
  /// \warning Overwrites nodes in ["fromNIdx-steps","fromNIdx") !
  /// \complexity O(1)
  inline void shift_nodes_up_(const NIdx fromNIdx, const NIdx toNIdx,
                              const Ind steps) noexcept {
    TRACE_IN((fromNIdx)(toNIdx)(steps));

    ASSERT(steps > 0, "Shifting nodes zero steps does nothing!");
    assert_valid_node_range(fromNIdx, toNIdx);
    assert_in_node_range(fromNIdx);
    assert_in_node_range(toNIdx-1);

    /// Copy nodes up in forward order to avoid overwritting:
    /// from first_node(toCIdx+1) to last_node(lastCell-1):
    for (NIdx nIdx = fromNIdx, endNIdx = toNIdx; nIdx != endNIdx; ++nIdx) {
      c()->copy_node_variables(nIdx, nIdx - steps);
    }

    TRACE_OUT();
  }

  /// \brief Shifts all nodes up \p steps times starting at the first node
  /// of the cell \p cIdx
  ///
  /// by updating the node indices for the cells in the range [cIdx,size), and
  /// shrinks the container size by \p steps nodes.
  ///
  /// \complexity O(M) where M is the #of _nodes_.
  ///
  inline void shift_cell_nodes_up_(const CIdx cIdx, const CIdx steps) noexcept {
    TRACE_IN((cIdx)(steps));
    ASSERT(steps > 0, "Shifting nodes zero steps does nothing!");
    assert_valid_node_range(first_node(cIdx), last_node(cIdx));
    ASSERT(node_size(cIdx) > 0, "Cell has no nodes!");
    ASSERT(cIdx > 0, "Cannot shift nodes up the zero-th cell!");

    /// Shift up nodes in range [first,last):
    const auto lastNode = node_size();
    const auto firstCellNode = first_node(cIdx);
    shift_nodes_up_(firstCellNode, lastNode, steps);

    /// Update cell node ranges: from (cIdx,0) ... (size(),1)
    // Dear optimizer from the future, we dont want to update the last
    // node (the invisible special one that lies inside nodes), so
    // we can't do array.tail(size()-cIdx).
    nodes_().array().segment(cIdx, size() - cIdx) -= steps;

    pop_node(steps);  ///< Shrink container.
    TRACE_OUT();
  }

  /// \brief Shifts node range [\p fromNIdx, \p toNIdx) down \p
  /// steps times
  ///
  /// \warning overwrites nodes in range (\p toNIdx, \p toNIdx+steps ] ! )
  ///
  /// \complexity O(1)
  ///
  inline void shift_nodes_down_(const NIdx fromNIdx, const NIdx toNIdx,
                                const Ind steps) noexcept {
    TRACE_IN((fromNIdx)(steps));

    ASSERT(steps > 0, "Shifting nodes zero steps does nothing!");
    ASSERT(toNIdx + steps < node_capacity(), "Out of memory!");
    assert_valid_node_range(fromNIdx, toNIdx);

    /// Copy nodes in reverse order to avoid overwritting
    // Note: fromNIdx can be zero, i.e. need to loop until nIdx >= endNIdx,
    // however --(nIdx=0) wraps around: nIdx = std::num_lim<Ind>::max() > 0!
    for (NIdx nIdx = toNIdx - 1, endNIdx = fromNIdx;
         nIdx >= endNIdx && nIdx != std::numeric_limits<NIdx>::max();
         --nIdx) {
      c()->copy_node_variables(nIdx, nIdx + steps);
    }

    TRACE_OUT();
  }

  /// \brief Shift all nodes down \p steps times starting at the first node of
  /// the cell \p cIdx
  ///
  /// by updating the node indices for the cells in the range
  /// [ \p cIdx, \p size), and grows the node container.
  ///
  /// \complexity O(M) where M is the #of _nodes_.
  ///
  inline void shift_cell_nodes_down_(const CIdx cIdx,
                                     const Ind steps) noexcept {
    TRACE_IN((cIdx)(steps));

    ASSERT(steps > 0, "Shifting nodes zero steps does nothing!");
    assert_valid_node_range(first_node(cIdx), last_node(cIdx));
    ASSERT(node_size(cIdx) > 0, "Cell has no nodes!");
    ASSERT(node_size() + steps < node_capacity(), "Out of memory!");

    /// Shift down nodes in range [first,last):
    const auto firstCellNode = first_node(cIdx);
    const auto lastNode = node_size();
    push_node(steps);  ///< Grows container.
    shift_nodes_down_(firstCellNode, lastNode, steps);

    /// Update cell node ranges: from (cIdx,0) ... (size(),1)
    // Dear optimizer from the future, we dont want to update the last
    // node (the invisible special one that lies inside nodes), so
    // we can't do array.tail(size()-cIdx).
    nodes_().array().segment(cIdx, size() - cIdx) += steps;
    TRACE_OUT();
  }

  /// \brief Resizes the number of nodes in a cell
  ///
  /// \warning after resizing, the value of the nodal variables is
  /// undefined. The value of the already existing nodes is presever _only if_
  /// the #of cell node grows or stays the same.
  ///
  /// \complexity O(M) where M is the #of _nodes_.
  ///
  inline void resize_cell_nodes_(const CIdx cIdx,
                                 const node_size_type newNoNodes) noexcept {
    TRACE_IN((cIdx)(newNoNodes));
    const auto noNodes = node_size(cIdx);
    const auto noNodesDifference = math::absdiff(noNodes, newNoNodes);

    if (noNodes < newNoNodes) {  ///< Grow cell
      shift_cell_nodes_down_(cIdx + 1, noNodesDifference);
      ASSERT(last_node(size()-1) == node_size(),
             "last's cell node != last node!");
      TRACE_OUT();
      return;
    } else if (noNodes > newNoNodes) {  ///< Shrink cell
      shift_cell_nodes_up_(cIdx + 1, noNodesDifference);
      ASSERT(last_node(size()-1) == node_size(),
             "last's cell node != last node!");
      TRACE_OUT();
      return;
    } else {  ///< noNodes == newNoNodes -> Resizing not necessary!
      ASSERT(last_node(size()-1) == node_size(),
             "last's cell node != last node!");
      TRACE_OUT();
      return;
    }

    ASSERT(false, "You should never get here!");
  }
  ///@}

  /// \brief Copies cells with fixed numbers of nodes
  ///
  /// \complexity O(1)
  ///
  void copy_cell_(const CIdx fromCIdx, const CIdx toCIdx,
                  tag::fixed_nodes) noexcept {
    TRACE_IN((fromCIdx)(toCIdx));
    c()->copy_cell_variables(fromCIdx, toCIdx);
    TRACE_OUT();
  }

  /// \brief Copies cells with variable numbers of nodes
  /// \warning Should be pretty slow.
  ///
  /// \complexity O(M) where M is the #of _nodes_!
  ///
  void copy_cell_(const CIdx fromCIdx, const CIdx toCIdx,
                  tag::variable_nodes) noexcept {
    TRACE_IN((fromCIdx)(toCIdx));
    ASSERT(fromCIdx != toCIdx, "Self copy not implemented!");

    /// Resize first, s.t. cell preserves its original state if resize fails!
    resize_cell_nodes_(toCIdx, node_size(fromCIdx));
    c()->copy_cell_variables(fromCIdx, toCIdx);
    copy_cell_nodes_(fromCIdx, toCIdx);
    TRACE_OUT();
  }

  /// \brief Copies all node variables from cell \p fromCIdx into \p toCIdx.
  ///
  /// \warning \p toCIdx is required to have exactly as many nodes as
  /// \p fromCIdx.
  ///
  /// \complexity O(1)
  ///
  inline void copy_cell_nodes_(const CIdx fromCIdx,
                               const CIdx toCIdx) noexcept {
    TRACE_IN((fromCIdx)(toCIdx));
    ASSERT(node_size(fromCIdx) == node_size(toCIdx),
           "Node range sizes are not equal!");

    const auto noCellNodes = node_size(fromCIdx);
    for (node_size_type i = 0; i < noCellNodes; ++i) {
      c()->copy_node_variables(first_node(fromCIdx) + i,
                               first_node(toCIdx) + i);
    }

    TRACE_OUT();
  }

  /// \brief Reset variables of all cells in range [fromCIdx,toCIdx)
  inline void reset_cells_(CIdx fromCIdx, const CIdx toCIdx) noexcept {
    for (; fromCIdx < toCIdx; ++fromCIdx) {
      c()->reset_cell(fromCIdx);
    }
  }
};


}  // namespace sequential

template<class T> using Sequential = sequential::Implementation<T>;

////////////////////////////////////////////////////////////////////////////////
}  // namespace container
}  // namespace hom3
////////////////////////////////////////////////////////////////////////////////
#undef assert_variable_node_container
#undef assert_in_cell_range
#undef assert_in_node_range
#undef assert_in_cell_node_range
#undef assert_valid_node_range
#undef ENABLE_DBG_
////////////////////////////////////////////////////////////////////////////////
#endif
