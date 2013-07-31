#ifndef HOM3_SOLVERS_FV_CONTAINER_HPP_
#define HOM3_SOLVERS_FV_CONTAINER_HPP_
/// Includes:
#include "../../containers/sequential.hpp"
/// Options:
#define ENABLE_DBG_ 0
#include "../../misc/dbg.hpp"
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
};
}} // container::sequential namespace

namespace solver { namespace fv {

/// \brief Extends container with global id functions
template<class CellContainer> struct GlobalIdFunctionality {
  inline const CellContainer* gidf_base() const {
    return static_cast<const CellContainer*>(this);
  }

  bool is_global(const Ind cellId) const {
    return is_valid<Ind>(gidf_base()->globalId(cellId));
  }

  inline RangeFilter<Ind> global() const {
    return RangeFilter<Ind>{[&](const Ind nId){ return is_global(nId); }};
  }

  /// \brief Returns [FilteredRange] of all global cell Ids.
  inline auto global_cells() const -> FRange<Ind> {
    return gidf_base()->all_cells() | global();
  }
};

template<SInd nd, SInd nvars>
struct Container : container::Sequential<Container<nd,nvars>>,
                   GlobalIdFunctionality<Container<nd,nvars>> {

  /// Aliases:
  friend container::Sequential< Container<nd,nvars> >;
  template<template <SInd> class V_, SInd nd_ = 1>
  using M = container::Matrix<Container,container::matrix::tag::Cell,V_,nd_>;

  /// Data:
  M<NumM,nvars> vars;
  M<NumM,nvars> rhs;
  M<IndM,2*nd>  nghbrs;
  M<NumM,nd>    xc;
  M<NumM,2*nd>  dx;
  M<NumM>       length;
  M<SIndM>      bcId;
  M<IndM>       globalId;

  /// Construction:
  Container(const Ind n) : container::Sequential<Container<nd,nvars>>(n,"fv_container"),
      vars(this,"vars"), rhs(this,"rhs"), nghbrs(this,"nghbrs"), xc(this,"xc"), dx(this,"xc"),
      length(this,"length"), bcId(this,"bcId") , globalId(this,"globalId") {
    std::cerr << "n cont: " << n << " capacity: " << this->capacity() << " size: " << this->size() << "\n";
  }

  inline void reset_cell(const Ind cId) {

    bcId(cId) = iSInd();
    globalId(cId) = iInd();
    length(cId) = iNum();
    for(SInd d = 0; d < nd; ++d) {
      xc(cId,d) = iNum();
    }
    for(SInd d = 0; d < 2*nd; ++d) {
      nghbrs(cId,d) = iInd();
      dx(cId,d) = iNum();
    }
    for(SInd v = 0; v < nvars; ++v) {
      rhs(cId,v) = iNum();
      vars(cId,v) = iNum();
    }
  }

 private:

  inline void copy_cell_variables(const Ind fromId, const Ind toId) {
    bcId(toId) = bcId(fromId);
    globalId(toId) = globalId(fromId);
    length(toId) = length(fromId);
    for(SInd d = 0; d < nd; ++d) {
      xc(toId,d) = xc(fromId,d);
    }
    for(SInd d = 0; d < 2*nd; ++d) {
      nghbrs(toId,d) = nghbrs(fromId,d);
      dx(toId,d) = dx(fromId,d);
    }
    for(SInd v = 0; v < nvars; ++v) {
      vars(toId,v) = vars(fromId,v);
      rhs(toId,v) = rhs(fromId,v);
    }
  }

};

/// Boilerplate: Value type
template<SInd nd, SInd nvars>
struct Value : container::sequential::ValueFacade<Container<nd,nvars>> {
  inline NumA<nvars>& vars()           { return vars_(); }
  inline Num&         vars(SInd d)     { return vars_(d); }
  inline NumA<nvars>& rhs()            { return rhs_(); }
  inline Num&         rhs(SInd d)      { return rhs_(d); }
  inline IndA<2*nd>&  nghbrs()         { return nghbrs_; }
  inline Ind&         nghbrs(SInd p)   { return nghbrs_(p); }
  inline NumA<2*nd>&  dx()             { return dx_; }
  inline Num&         dx(SInd p)       { return dx_(p); }
  inline NumA<nd>&    xc()             { return xc_; }
  inline Num&         xc(SInd d)       { return xc_(d); }
  inline SInd&        bcId()           { return bcId_; }
  inline Ind&         globalId()       { return globalId_; }
  inline Num&         length()         { return length_; }

  template<class Value1, class Value2>
  static inline void swap_values(Value1& lhs, Value2& rhs) {
    std::swap(lhs.bcId(),rhs.bcId());
    std::swap(lhs.globalId(),rhs.globalId());
    std::swap(lhs.length(),rhs.length());
    for(SInd p = 0; p < 2*nd; ++p) {
      std::swap(lhs.nghbrs(p),rhs.nghbrs(p));
      std::swap(lhs.dx(p),rhs.dx(p));
    }
    for(SInd d = 0; d < nd; ++d) {
      std::swap(lhs.xc(d),rhs.xc(d));
    }
    for(SInd v = 0; v < nvars; ++v) {
      std::swap(lhs.vars(v),rhs.vars(v));
      std::swap(lhs.rhs(v),rhs.rhs(v));
    }
  }

  template<class Value1, class Value2>
  static inline void copy_values(Value1& lhs, Value2& rhs) {
    lhs.bcId() = rhs.bcId();
    lhs.globalId() = rhs.globalId();
    lhs.length() = rhs.length();
    for(SInd p = 0; p < 2*nd; ++p) {
      lhs.nghbrs(p) = rhs.nghbrs(p);
      lhs.dx(p) = rhs.dx(p);
    }
    for(SInd d = 0; d < nd; ++d) {
      lhs.xc(d) = rhs.xc(d);
    }
    for(SInd v = 0; v < nvars; ++v) {
      lhs.vars(v),rhs.vars(v);
      lhs.rhs(v),rhs.rhs(v);
    }

  }
 private:
  NumA<nvars> vars_, rhs_;
  SInd bcId_;
  Ind globalId_;
  IndA<2*nd> nghbrs_;
  NumA<2*nd> dx_;
  NumA<nd> xc_;
  Num length_;
};

/// Boilerplate: Reference Type
template<SInd nd, SInd nvars>
struct Reference : container::sequential::ReferenceFacade<Container<nd,nvars>> {
  using Base = container::sequential::ReferenceFacade<Container<nd,nvars>>;
  using Base::c;
  using Base::index;
  using Base::operator=;
  Reference<nd,nvars>& operator=(Reference rhs_) { return Base::operator=(rhs_); }

  inline Num&  vars(SInd d = 0) { return c()->vars(index(),d); }
  inline Num&  rhs(SInd d = 0) { return c()->rhs(index(),d); }
  inline Ind&  nghbrs(SInd d) { return c()->nghbrs(index(),d); }
  inline Num&  dx(SInd d) { return c()->dx(index(),d); }
  inline Num&  xc(SInd d) { return c()->xc(index(),d); }
  inline SInd& bcId() { return c()->bcId(index()); }
  inline Ind&  globalId() { return c()->globalId(index()); }
  inline Num&  length() { return c()->length(index()); }

  Reference() : Base()  {}
  Reference(Container<nd,nvars>* c, Ind i) : Base(c,i) {}
};

}} // solver::fv namespace

#undef ENABLE_DBG_
////////////////////////////////////////////////////////////////////////////////
#endif
