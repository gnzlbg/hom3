#ifndef HOM3_SOLVERS_FV_CONTAINER_HPP_
#define HOM3_SOLVERS_FV_CONTAINER_HPP_
/// Includes:
#include <algorithm>
#include "../../containers/sequential.hpp"
/// Options:
#define ENABLE_DBG_ 0
#include "../../misc/dbg.hpp"
////////////////////////////////////////////////////////////////////////////////
namespace hom3 {
////////////////////////////////////////////////////////////////////////////////

/// Forward declarations:
namespace solver { namespace fv {

template<SInd nd,SInd nvars> struct Container;
template<SInd nd,SInd nvars> struct Reference;
template<SInd nd,SInd nvars> struct Value;

}} // solver::fv namespace

/// Cell traits specialization:
namespace container { namespace sequential {
template<SInd nd,SInd nvars> struct traits<solver::fv::Container<nd,nvars>> {
  using container_type    = tag::fixed_nodes;
  using value_type        = typename solver::fv::Value<nd,nvars>;
  using reference         = typename solver::fv::Reference<nd,nvars>;
  using cell_index_type   = CellIdx;
};
}} // container::sequential namespace

namespace solver { namespace fv {

/// \brief Extends container with node id functions
template<class CellContainer> struct NodeIdsFunctionality {
  inline const CellContainer* gidf_base() const {
    return static_cast<const CellContainer*>(this);
  }

  bool is_node(const Ind cellId) const {
    return is_valid<Ind>(gidf_base()->node_idx(cellId));
  }

  inline RangeFilter<Ind> is_node() const {
    return RangeFilter<Ind>{[&](const Ind nId){ return is_node(nId); }};
  }

  /// \brief Returns [FilteredRange] of all global cell Ids.
  inline auto global_cells() const -> FRange<Ind> {
    return gidf_base()->all_cells() | is_node();
  }
};

template<SInd nd, SInd nvars>
struct Container : container::Sequential<Container<nd,nvars>>,
                   NodeIdsFunctionality<Container<nd,nvars>> {

  /// Aliases:
  friend container::Sequential< Container<nd,nvars> >;
  friend struct Reference<nd,nvars>;
  template<template <SInd> class V_, SInd nd_ = 1>
  using M = container::Matrix<Container,container::matrix::tag::Cell,V_,CellIdx,SInd,nd_>;

  /// Construction:
  Container(const Ind n) : container::Sequential<Container<nd,nvars>>(n,"fv_container"),
      lhs(this,"lhs"), rhs(this,"rhs"), neighbors(this,"neighbors"),
      x_center(this,"x"), distances(this,"distances"), length(this,"length"),
      bc_idx(this,"bc_idx") , node_idx(this,"node_idx") {
    // std::cerr << "n cont: " << n << " capacity: " << this->capacity()
    //           << " size: " << this->size() << "\n";
  }

  /// Data:
  M<NumM,nvars>    lhs;
  M<NumM,nvars>    rhs;
  M<CellIdxM,2*nd> neighbors;
  M<NumM,nd>       x_center;
  M<NumM,2*nd>     distances;
  M<NumM>          length;
  M<SIndM>         bc_idx;
  M<NodeIdxM>      node_idx;

  inline void reset_cell(const CellIdx cId) {
    bc_idx(cId) = invalid<SInd>();
    node_idx(cId) = invalid<NodeIdx>();
    length(cId) = invalid<Num>();
    for(SInd d = 0; d < nd; ++d) {
      x_center(cId,d) = invalid<Num>();
    }
    for(SInd d = 0; d < 2*nd; ++d) {
      neighbors(cId,d) = invalid<CellIdx>();
      distances(cId,d) = invalid<Num>();
    }
    for(SInd v = 0; v < nvars; ++v) {
      rhs(cId,v) = invalid<Num>();
      lhs(cId,v) = invalid<Num>();
    }
  }

  inline void copy_cell_variables(const CellIdx fromId, const CellIdx toId) {
    bc_idx(toId) = bc_idx(fromId);
    node_idx(toId) = node_idx(fromId);
    length(toId) = length(fromId);
    for(SInd d = 0; d < nd; ++d) {
      x_center(toId,d) = x_center(fromId,d);
    }
    for(SInd d = 0; d < 2*nd; ++d) {
      neighbors(toId,d) = neighbors(fromId,d);
      distances(toId,d) = distances(fromId,d);
    }
    for(SInd v = 0; v < nvars; ++v) {
      lhs(toId,v) = lhs(fromId,v);
      rhs(toId,v) = rhs(fromId,v);
    }
  }

};

/// Boilerplate: Value type
template<SInd nd, SInd nvars>
struct Value : container::sequential::ValueFacade<Container<nd,nvars>> {
  inline NumA<nvars>& lhs()             { return lhs_;          }
  inline Num&         lhs(SInd d)       { return lhs_(d);       }
  inline NumA<nvars>& rhs()             { return rhs_;          }
  inline Num&         rhs(SInd d)       { return rhs_(d);       }
  inline CellIdxA<2*nd>&  neighbors()   { return neighbors_;    }
  inline CellIdx&     neighbors(SInd p) { return neighbors_(p); }
  inline NumA<2*nd>&  distances()       { return distances_;    }
  inline Num&         distances(SInd p) { return distances_(p); }
  inline NumA<nd>&    x()               { return x_;            }
  inline Num&         x(SInd d)         { return x_(d);         }
  inline SInd&        bc_idx()          { return bc_idx_;       }
  inline NodeIdx&     node_idx()        { return node_idx_;     }
  inline Num&         length()          { return length_;       }

  template<class Value1, class Value2>
  static inline void swap_values(Value1& a, Value2& b) {
    using std::swap;
    swap(a.bc_idx(),b.bc_idx());
    swap(a.node_idx(),b.node_idx());
    swap(a.length(),b.length());
    for(SInd p = 0; p < 2*nd; ++p) {
      swap(a.neighbors(p),b.neighbors(p));
      swap(a.distances(p),b.distances(p));
    }
    for(SInd d = 0; d < nd; ++d) {
      swap(a.x(d),b.x(d));
    }
    for(SInd v = 0; v < nvars; ++v) {
      swap(a.lhs(v),b.lhs(v));
      swap(a.rhs(v),b.rhs(v));
    }
  }

  template<class Value1, class Value2>
  static inline void copy_values(Value1& a, Value2& b) {
    a.bc_idx() = b.bc_idx();
    a.node_idx() = b.node_idx();
    a.length() = b.length();
    for(SInd p = 0; p < 2*nd; ++p) {
      a.neighbors(p) = b.neighbors(p);
      a.distances(p) = b.distances(p);
    }
    for(SInd d = 0; d < nd; ++d) {
      a.x(d) = b.x(d);
    }
    for(SInd v = 0; v < nvars; ++v) {
      a.lhs(v),b.lhs(v);
      a.rhs(v),b.rhs(v);
    }

  }
 private:
  NumA<nvars> lhs_, rhs_;
  SInd bc_idx_;
  NodeIdx node_idx_;
  CellIdxA<2*nd> neighbors_;
  NumA<2*nd> distances_;
  NumA<nd> x_;
  Num length_;
};

/// Boilerplate: Reference Type
template<SInd nd, SInd nvars>
struct Reference : container::sequential::ReferenceFacade<Container<nd,nvars>> {
  using Base = container::sequential::ReferenceFacade<Container<nd,nvars>>;
  using Base::c;
  using Base::index;
  using Base::operator=;
  using Base::Base;
  Reference<nd,nvars>& operator=(Reference rhs_) { return Base::operator=(rhs_); }

  inline Num& lhs(SInd d = 0)       { return c()->lhs(index(),d);       }
  inline Num& rhs(SInd d = 0)       { return c()->rhs(index(),d);       }
  inline CellIdx& neighbors(SInd d) { return c()->neighbors(index(),d); }
  inline Num& distances(SInd d)     { return c()->distances(index(),d); }
  inline Num& x(SInd d)             { return c()->x_center(index(),d);  }
  inline SInd& bc_idx()             { return c()->bc_idx(index());      }
  inline NodeIdx& node_idx()        { return c()->node_idx(index());    }
  inline Num& length()              { return c()->length(index());      }
};

}} // solver::fv namespace

////////////////////////////////////////////////////////////////////////////////
} // hom3 namespace
////////////////////////////////////////////////////////////////////////////////
#undef ENABLE_DBG_
////////////////////////////////////////////////////////////////////////////////
#endif
