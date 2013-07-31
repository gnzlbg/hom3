#ifndef HOM3_CONTAINER_SEQUENTIAL_IMPLEMENTATION_HPP_
#define HOM3_CONTAINER_SEQUENTIAL_IMPLEMENTATION_HPP_
////////////////////////////////////////////////////////////////////////////////
/// \file \brief Implementation of the sequential container
////////////////////////////////////////////////////////////////////////////////
/// Includes:
#include <type_traits>
/// Options:
#define ENABLE_DBG_ 0
#include "../../misc/dbg.hpp"
////////////////////////////////////////////////////////////////////////////////

#define assert_variable_node_container()                                \
  static_assert(std::is_same<container_type,tag::variable_nodes>::value, \
                "Function defined only for variable_node_size containers!")

namespace container { namespace sequential {

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
  using size_type      = Ind;

  /// \name Constructors
  ///@{
  Implementation(const Ind n, const std::string name) : Implementation(n,0,name) {}
  Implementation(const Ind ne, const Ind nn, const std::string name)
      : maxCellSize_(ne), maxNodeSize_(nn), cellSize_(0),
        nodeSize_(0), nodes_(c(),"nodes"), name_(name) {
    TRACE_IN((ne)(nn));

    if (capacity() == 0) {
      TERMINATE("Empty container!");
    }

    if (std::is_same<container_type,tag::variable_nodes>::value
        && node_capacity() < capacity()) {
      TERMINATE("Less nodes than cells!");
    }

    TRACE_OUT();
  }

  Implementation(const Implementation& other)
      : maxCellSize_(other.capacity()), maxNodeSize_(other.node_capacity()),
        cellSize_(other.size()), nodeSize_(other.node_size()), nodes_(other.nodes_),
        name_(other.name_) {
    TRACE_IN_();

    init_nodes_(container_type());

    TRACE_OUT();
  }

  Implementation(Implementation&& other)
      : maxCellSize_(other.capacity()), maxNodeSize_(other.node_capacity()),
        cellSize_(other.size()), nodeSize_(other.node_size()), nodes_(other.nodes_)
  {
    TRACE_IN_();
    TRACE_OUT();
  }


  /// \brief Initializes a single column type
  template<class Variable> void initialize(Variable&& variable) {
    variable.init(c());
  }
  /// \brief Initialize the column types by calling init for each member  (each
  /// member's init function takes the container pointer and uses it to check
  /// ownership)
  template<class Variable, class... Variables>
  void initialize (Variable&& variable, Variables&&... variables) {
    initialize(std::forward<Variable>(variable));
    initialize(std::forward<Variables>(variables)...);
  }
  ///@}

  /// \name Size functions
  ///@{

  /// \brief Returns the maximum #of cells that can be stored in the container
  inline Ind  capacity()      const { return maxCellSize_;   }
  /// \brief Returns the maximum #of nodes that can be stored in the container
  inline Ind  node_capacity() const { return maxNodeSize_;   }
  /// \brief Is the container empty (i.e. does it have no cells?)
  inline bool empty()         const { return cellSize_ == 0; }
  /// \brief Returns the current #of cells stored in the container
  inline Ind  size()          const { return cellSize_;      }
  /// \brief Returns the current #of nodes stored in the container
  inline Ind  node_size()     const
  { assert_variable_node_container(); return nodeSize_;      }
  /// \brief Returns the current #of nodes of the element \p cellId
  inline Int node_size(const Ind cellId) { ///< Returns #nodes of "cellId"
    assert_variable_node_container();
    ASSERT(first_node(cellId) <= last_node(cellId),"Invalid cell node range!");
    return last_node(cellId) - first_node(cellId);
  }
  ///@}

  /// \name Cell iterators
  ///@{

  /// \brief Returns a cell iterator to the first container element
  inline iterator        begin()       { return {c(),0};      }
  /// \brief Returns a const cell iterator to the first container element
  inline const_iterator cbegin() const { return {c(),0};      }

  /// \brief Returns a cell iterator to the one-past-the-end container element
  inline iterator          end()       { return {c(),size()}; }
  /// \brief Returns a const cell iterator to the one-past-the-end container element
  inline const_iterator   cend() const { return {c(),size()}; }

  ///@}

  /// \name Indexed cell access
  ///@{

  /// \brief Returns index of the first container element (always 0)
  inline Ind front() const { return 0; }
  /// \brief Returns index of the last container element
  inline Ind back()  const { return size() - 1; }
  /// \brief Random access: returns a reference to the element \p i
  inline reference operator[](const Ind i) {
    ASSERT(i < size(),"Index is out of bounds.");
    return {c(),i};
  }
  inline reference operator[](const Ind i) const {
    ASSERT(i < size(),"Index is out of bounds.");
    return {c(),i};
  }
  ///@}

  /// \name Indexed ranges
  ///@{

  /// \brief Returns [IndRange] of _all_ container cell Ids.
  inline auto all_cells() const -> Range<Ind> { return Range<Ind>{Ind(0),size()}; }
  ///@}

  ///@}

  /// \name Indexed node access
  ///@{

  /// \brief Returns index of the first node of the cell \p cellId
  inline Ind first_node(Ind cellId = 0)  { return nodes_(cellId);        }

  /// \brief Returns index of the last node of the cell \p cellId
  ///
  /// Note: In the last cell "last_node(size()-1)" calls "nodes_(size())".
  /// Although the cell with cellId == size() does not exists, the container has
  /// an extra _magic_ node that marks the end of the last cells' nodes.
  ///
  ///@{
  inline Ind last_node (  Ind cellId  )  { return nodes_(cellId + 1);    }
  inline Ind last_node (              )  { return last_node(size() - 1); }
  ///@}
  ///@}

  /// \name Append/delete cells

  ///@{

  /// \brief Inserts \p noCells with \p noNodesPerCell at the container end
  /// and resets their variables
  void push_cell(const Ind noCells = 1, const Int noNodesPerCell = 1) {
    push_cell_(noCells,noNodesPerCell,container_type());
    reset_cells_(size() - noCells, size());
  }
  /// \brief Removes last \p i cells from the container
  void pop_cell(const Ind i = 1) {
    pop_cell_(i,container_type());
  }
  ///@}

  /// \name Append/delete nodes
  ///@{

  /// \brief Inserts \p i nodes at the end of the container
  void push_node(const Ind i = 1) {
    node_size_() += i;
    last_node_() += i;
    ASSERT(node_size() <= node_capacity(), "Container out of memory.");
  }
  /// \brief Removes \p i nodes at the end of the container
  void pop_node(const Ind i = 1) {
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
  Ind erase_remove_if(Predicate&& p) {
    const Ind last = size();
    Ind first = algorithm::find_if(0,last,p);

    if(first == last) { return first; }

    Ind next = first; ++first;
    for(; first != last; ++first)
      if(!p(first)) {
        copy_cell(first,next);
        ++next;
      }

    pop_cell(last-next);
    return next;
  }

  /// \brief Copies the element at \p fromCellId to position \p toCellId
  ///
  /// \complexity O(1) for fixed_nodecontainers
  /// \complexity O(N) for variable_nodes containers
  void copy_cell(const Ind fromCellId, const Ind toCellId) {
    ASSERT(fromCellId != toCellId, "Suspicious | Untested behaviour");
    copy_cell_(fromCellId,toCellId,container_type());
  }

  ///@}

  /// \name Extra functionality
  ///@{

  /// \brief Swap containers
  /// \todo should swap be limited to containers of equal maxCell/Node sizes?
  friend void swap(Cells& first, Cells& second) {
    TRACE_IN_();

    /// Swaps container variables:
    std::swap(first.maxCellSize_ , second.maxCellSize_ );
    std::swap(first.maxNodeSize_ , second.maxNodeSize_ );
    std::swap(first.cellSize_    , second.cellSize_    );
    std::swap(first.nodeSize_    , second.nodeSize_    );
    std::swap(first.nodes_       , second.nodes_       );

    /// User defined function for swapping user-defined variables:
    swap_containers(first,second);

    TRACE_OUT();
  }
  ///@}

  std::string name() const { return name_; }

 private:

  // Note: max sizes aren't const to allow swapping containers of different sizes
  // \todo remove swapping of containers and make them const
  Ind maxCellSize_; ///< max #of cells + 1
  Ind maxNodeSize_; ///< max #of nodes + 1

  Ind cellSize_; ///< #of cells + 1
  Ind nodeSize_; ///< #of nodes + 1

  /// Container of node indicies
  NodeIndices<Cells> nodes_;

  const std::string name_;

  /// \brief Returns a pointer to the underlying container (CRTP)
  inline Cells* c() { return static_cast<Cells*>(this); }

  /// \brief Generic initialization of both fixed and variable node containers
  ///@{
  void init_nodes_(tag::fixed_nodes ) {}
  void init_nodes_(tag::variable_nodes ) { initialize(nodes_); }
  ///@}

  /// \name Allow to modify the container private variables
  ///@{
  inline Ind& size_()      { return cellSize_; }
  inline Ind& node_size_() { return nodeSize_; }
  inline Ind& first_node_(Ind cellId = 0) { return nodes_(cellId);         }
  inline Ind& last_node_ (Ind cellId    ) { return nodes_(cellId + 1);     }
  inline Ind& last_node_ (              ) { return last_node_(size() - 1); }
  ///@}

  /// Finds the cellId of a given node (slow: for debugging purposes only!)
  inline Ind findCellWithNode(Ind nodeId) {
    for(Ind cId = 0; cId < size(); ++cId)
      if(first_node(cId) <= nodeId && last_node(cId) > nodeId)
        return cId;
    return size();
  }

  /// \name Implementation details of append/delete functions
  ///@{
  void push_cell_(const Ind i, Ind, tag::fixed_nodes) {
    size_() += i;
    ASSERT(size() <= capacity(), "Container out of memory.");
  }

  void push_cell_(const Ind noCells, const Int noNodesPerCell, tag::variable_nodes) {
    Ind cellId = size(), nodeId = node_size();
    size_() += noCells;
    node_size_() += noCells * noNodesPerCell;
    for(; cellId < size(); ++cellId) {
      first_node_(cellId) = nodeId;
      nodeId += noNodesPerCell;
    }
    last_node_() = node_size();
    ASSERT(size() <= capacity(), "Container out of memory.");
  }

  void pop_cell_(const Ind i, tag::fixed_nodes){
    ASSERT(!empty(), "Container is already empty!");
    ASSERT(size() - i > 0, "Container doesn't have enough cells to pop!");
    size_() -= i;
  }

  void pop_cell_(const Ind i, tag::variable_nodes){
    ASSERT(!empty(), "Container is already empty!");
    ASSERT(size() - i > 0, "Container doesn't have enough cells to pop!");
    last_node_() = last_node(size() - i - 1);
    node_size_() = last_node();
    size_() -= i;
  }

  /// \brief Shifts node range ["fromNodeId","toNodeId") up "steps" times
  ///
  /// \warning Overwrites nodes in ["fromNodeId-steps","fromNodeId") !
  /// \complexity O(1)
  inline void shift_nodes_up_(const Ind fromNodeId, const Ind toNodeId, const Ind steps) {
    TRACE_IN((fromNodeId)(toNodeId)(steps));

    ASSERT(steps > 0, "Shifting nodes zero steps does nothing!");
    ASSERT(fromNodeId < toNodeId, "Node out of range!");
    ASSERT(fromNodeId != toNodeId, "Shifting empty range!");
    ASSERT(fromNodeId <= node_size(), "Node out of range!");
    ASSERT(toNodeId <= node_size(), "Node out of range!");

    /// Copy nodes up in forward order to avoid overwritting:
    /// from first_node(toCellId+1) to last_node(lastCell-1):
    for(Ind nodeId = fromNodeId,
         endNodeId = toNodeId; nodeId != endNodeId; ++nodeId)
      c()->copy_node_variables(nodeId,nodeId-steps);

    TRACE_OUT();
  }

  /// \brief Shifts all nodes up \p steps times starting at the first node
  /// of the cell \p cellId
  ///
  /// by updating the node indices for the cells in the range [cellId,size), and
  /// shrinks the container size by \p steps nodes.
  ///
  /// \complexity O(M) where M is the #of _nodes_.
  ///
  inline void shift_cell_nodes_up_(const Ind cellId, const Ind steps) {
    TRACE_IN((cellId)(steps));

    ASSERT(steps > 0, "Shifting nodes zero steps does nothing!");
    ASSERT(first_node(cellId) < last_node(cellId),"Invalid cell node range!");
    ASSERT(node_size(cellId) > 0, "Cell has no nodes!");
    ASSERT(cellId > 0,"Cannot shift up nodes from the zeroth cell!");

    /// Shift up nodes in range [first,last):
    const Ind lastNode = node_size();
    const Ind firstCellNode = first_node(cellId);

    shift_nodes_up_(firstCellNode, lastNode, steps);

    /// Update cell node ranges: from (cellId,0) ... (size(),1)
    // Dear optimizer from the future, we dont want to update the last
    // node (the invisible special one that lies inside nodes), so
    // we can't do array.tail(size()-cellId).
    nodes_().array().segment(cellId, size() - cellId) -= steps;

    pop_node(steps); ///< Shrink container.

    TRACE_OUT();
  }

  /// \brief Shifts node range [\p fromNodeId, \p toNodeId) down \p
  /// steps times
  ///
  /// \warning overwrites nodes in range (\p toNodeId, \p toNodeId+steps ] ! )
  ///
  /// \complexity O(1)
  ///
  inline void shift_nodes_down_(const Ind fromNodeId, const Ind toNodeId, const Ind steps) {
    TRACE_IN((fromNodeId)(steps));

    //DBGV((toNodeId + steps)(node_capacity()));
    ASSERT(steps > 0, "Shifting nodes zero steps does nothing!");
    ASSERT(toNodeId + steps < node_capacity(),"Out of memory!");
    ASSERT(fromNodeId < toNodeId, "Node out of range!");

    /// Copy nodes in reverse order to avoid overwritting
    // Note: fromNodeId can be zero, i.e. need to loop until nodeId >= endNodeId,
    // however --(nodeId=0) wraps around: nodeId = std::num_lim<Ind>::max() > 0!
    for(Ind nodeId = toNodeId - 1, endNodeId = fromNodeId;
        nodeId >= endNodeId && nodeId != std::numeric_limits<Ind>::max(); --nodeId)
      c()->copy_node_variables(nodeId, nodeId + steps);

    TRACE_OUT();
  }

  /// \brief Shift all nodes down \p steps times starting at the first node of the
  /// cell \p cellId
  ///
  /// by updating the node indices for the cells in the range
  /// [ \p cellId, \p size), and grows the node container.
  ///
  /// \complexity O(M) where M is the #of _nodes_.
  ///
  inline void shift_cell_nodes_down_(const Ind cellId, const Ind steps) {
    TRACE_IN((cellId)(steps));

    ASSERT(steps > 0, "Shifting nodes zero steps does nothing!");
    ASSERT(first_node(cellId) <= last_node(cellId),"Invalid cell node range!");
    ASSERT(node_size(cellId) > 0, "Cell has no nodes!");
    ASSERT(node_size() + steps < node_capacity(), "Out of memory!");

    /// Shift down nodes in range [first,last):
    const Ind firstCellNode = first_node(cellId);
    const Ind lastNode = node_size();

    push_node(steps); ///< Grows container.

    shift_nodes_down_(firstCellNode,lastNode,steps);

    /// Update cell node ranges: from (cellId,0) ... (size(),1)
    // Dear optimizer from the future, we dont want to update the last
    // node (the invisible special one that lies inside nodes), so
    // we can't do array.tail(size()-cellId).
    nodes_().array().segment(cellId, size() - cellId) += steps;
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
  inline void resize_cell_nodes_(const Ind cellId, const Ind newNoNodes) {
    TRACE_IN((cellId)(newNoNodes));

    const Ind noNodes = node_size(cellId);

    Ind noNodesDifference = math::absdiff(noNodes,newNoNodes);

    if (noNodes < newNoNodes) { ///< Grow cell

      shift_cell_nodes_down_(cellId+1,noNodesDifference);

      ASSERT(last_node(size()-1) == node_size(),"last's cell node != last node!");
      TRACE_OUT();
      return;

    } else if (noNodes > newNoNodes) { ///< Shrink cell

      shift_cell_nodes_up_(cellId+1,noNodesDifference);

      ASSERT(last_node(size()-1) == node_size(),"last's cell node != last node!");
      TRACE_OUT();
      return;

    } else { ///< noNodes == newNoNodes -> Resizing not necessary!
      ASSERT(last_node(size()-1) == node_size(),"last's cell node != last node!");
      TRACE_OUT();
      return;
    }

    ASSERT(false,"You should never get here!");
  }
  ///@}

  /// \brief Copies cells with fixed numbers of nodes
  ///
  /// \complexity O(1)
  ///
  void copy_cell_(const Ind fromCellId, const Ind toCellId, tag::fixed_nodes) {
    TRACE_IN((fromCellId)(toCellId));

    c()->copy_cell_variables(fromCellId,toCellId);

    TRACE_OUT();
  }

  /// \brief Copies cells with variable numbers of nodes
  /// \warning Should be pretty slow.
  ///
  /// \complexity O(M) where M is the #of _nodes_!
  ///
  void copy_cell_(const Ind fromCellId, const Ind toCellId, tag::variable_nodes) {
    TRACE_IN((fromCellId)(toCellId));

    ASSERT(fromCellId != toCellId, "Self copy not implemented!");

    /// Resize first, s.t. cell preserves its original state if resize fails!
    resize_cell_nodes_(toCellId,node_size(fromCellId));

    c()->copy_cell_variables(fromCellId,toCellId);

    copy_cell_nodes_(fromCellId,toCellId);

    TRACE_OUT();
  }

  /// \brief Copies all node variables from cell \p fromCellId into \p toCellId.
  ///
  /// \warning \p toCellId is required to have exactly as many nodes as
  /// \p fromCellId.
  ///
  /// \complexity O(1)
  ///
  inline void copy_cell_nodes_(const Ind fromCellId, const Ind toCellId) {
    TRACE_IN((fromCellId)(toCellId));

    ASSERT(node_size(fromCellId) == node_size(toCellId),      \
           "Node range sizes are not equal!");

    Ind noCellNodes = node_size(fromCellId);
    for(Ind i = 0; i < noCellNodes; ++i) {
      c()->copy_node_variables(first_node(fromCellId)+i,first_node(toCellId)+i);
    }

    TRACE_OUT();
  }

  /// \brief Reset variables of all cells in range [fromCellId,toCellId)
  inline void reset_cells_(Ind fromCellId, const Ind toCellId) {
    for(; fromCellId < toCellId; ++fromCellId) {
      c()->reset_cell(fromCellId);
    }
  }
};


} // namespace sequential

template<class T> using Sequential = sequential::Implementation<T>;

} // namespace container
#undef ENABLE_DBG_
////////////////////////////////////////////////////////////////////////////////
#endif
