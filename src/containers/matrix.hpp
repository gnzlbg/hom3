#ifndef HOM3_CONTAINER_MATRIX_HPP_
#define HOM3_CONTAINER_MATRIX_HPP_
////////////////////////////////////////////////////////////////////////////////
/// Includes:
#include "globals.hpp"
#include <type_traits>
#include <algorithm>
#include <string>
////////////////////////////////////////////////////////////////////////////////

namespace hom3 { namespace container {

/// \todo rename tags/move them somewhere else/replace tag dispatch with crtp.

namespace matrix {

namespace tag {
struct Cell {};
struct Node {};
struct NodeIndices {};
struct unknown {};
struct eigen {};
} // tag namespace

/// Column type traits
template<class T> struct traits {
  using value_type      = typename T::value_type;
  using reference       = typename T::reference;
  using const_reference = typename T::const_reference;
  using type            = tag::unknown;
};

template<class T, int r, int c, int t, int p1, int p2>
struct traits<Eigen::Matrix<T,r,c,t,p1,p2>> {
 using value_type      = typename Eigen::Matrix<T,r,c,t,p1,p2>::Scalar;
 using reference       = value_type      &;
 using const_reference = value_type const&;
 using type            = tag::eigen;
};

} // matrix namespace

/// \brief Stores multi-dimensional cell/nodal variables using the container \p C.
///
/// Wraps different types of containers (e.g. Eigen::Matrix, std::vector, ...)
/// providing the same interface
///
/// Template parameters:
/// C = Container Type
/// T = Variable semantics: cell/nodal variables
/// V<i> = Basic type: i columns of V type (variables_traits available)
/// nd_ = number of columns -> V<nd_>
template<class C, class T, template <SInd> class V, class RowIdx, class ColIdx, SInd nd_ = 1> struct Matrix {

  /// #of columns
  inline static constexpr SInd nd(){ return nd_; }

  using This             = Matrix;
  using container        = V<nd()>; ///< Underlying container type
  using container_trait  = matrix::traits<container>;
  using reference        = typename container_trait::reference;
  using const_reference  = typename container_trait::const_reference;
  using container_type   = typename container_trait::type;

  Matrix() : c_(nullptr), data_(0,nd()), name_("unknown") {}
  Matrix(const C* t, std::string name) : c_(t), data_(capacity_(),nd()), name_(name) {}
  Matrix(const This& other) : c_(nullptr), data_(other()), name_(other.name()) {}
  Matrix(This&& other) : Matrix() {
    using std::swap;
    swap(*this, other);
  }
  This& operator=(This other) {
    using std::swap;
    swap(*this, other);
    return *this;
  }

  /// \brief Access element at index (\p row_i, \p col_j)
  inline reference       operator()(const RowIdx row_i, const ColIdx col_j = 0) {
    ASSERT(preconditions_(row_i,col_j),"");
    return data_(primitive_cast(row_i),primitive_cast(col_j));
  }
  inline const_reference operator()(const RowIdx row_i, const ColIdx col_j = 0) const {
    ASSERT(preconditions_(row_i,col_j),"");
    return const_reference(data_(primitive_cast(row_i),primitive_cast(col_j)));
  }

  /// \brief Swaps data without swapping ownership:
  /// \todo swapping container ownership is easier and faster
  friend void swap(This& a, This& b) {
    using std::swap;
    swap(a.data_, b.data_);
    swap(a.name_,b.name_);
  }

  inline void init(const C* t) { c_ = t; data_.resize(capacity_(),nd()); }

  /// \brief Acces the underlying container type
  inline       container& operator()()       { return data_; }
  inline const container& operator()() const { return data_; }

  /// \brief Flips all bytes (work with containers of bools only!)
  inline void flip() { data_.flip(); } ///< flips all!
  inline void flip(const SInd d) { data_.flip(d); } //< flips all in dim nd_

  inline std::string name() const { return name_; }

 private:
  const C* c_;
  container data_;
  std::string name_;

  inline std::string path() const { return name() + " in container " + c_->name(); }

  inline bool preconditions_(const RowIdx i, const ColIdx d) const {
    ASSERT(data_.size() != 0,"Uninitialized variables: " << path() << " has size = 0!");
    ASSERT(primitive_cast(i) < size_(),
           "Index " << i << " is out of bounds: " << path() << " has size = " << size_() << ".");
    ASSERT(primitive_cast(i) < capacity_(),
           "Index " << i << " in " << path() << " has capacity = " << capacity_() << ".");
    ASSERT(d < nd(),"Column index " << d << " is out of bounds: " << path() << " has " << nd() << " columns.");
    return true;
  }

  // Note: If the list of cases below grows, maybe move behavior into a veneer?

  /// \brief Container capacity
  inline Ind capacity_() const { return capacity_impl_(T()); }
  inline Ind capacity_impl_(matrix::tag::Cell)        const { return c_->capacity();      }
  inline Ind capacity_impl_(matrix::tag::Node)        const { return c_->node_capacity(); }
  inline Ind capacity_impl_(matrix::tag::NodeIndices) const { return c_->capacity() + 1;  }

  /// \brief Container size
  inline Ind size_() const { return size_impl_(T()); }
  inline Ind size_impl_(matrix::tag::Cell)        const { return c_->size();      }
  inline Ind size_impl_(matrix::tag::Node)        const { return c_->node_size(); }
  inline Ind size_impl_(matrix::tag::NodeIndices) const { return c_->size() + 1;  }

 public:
  inline auto row(const RowIdx i) -> decltype(std::declval<This>().data_.row(primitive_cast(i))) {
    return data_.row(primitive_cast(i));
  }
  inline auto col(const ColIdx i) -> decltype(data_.col(primitive_cast(i))) {
    return data_.col(primitive_cast(i));
  }
  inline auto row(const RowIdx i) const -> decltype(std::declval<const This>().data_.row(primitive_cast(i))) {
    return data_.row(primitive_cast(i));
  }
  inline auto col(const ColIdx i) const -> decltype(data_.col(primitive_cast(i))) {
    return data_.col(primitive_cast(i));
  }

};

// Alias (sugar) for array of node indices
template<class C,class RowIdx,class ColIdx = SInd>
using NodeIndices = Matrix<C,matrix::tag::NodeIndices,IndM,RowIdx,ColIdx,1>;
////////////////////////////////////////////////////////////////////////////////
}} // hom3::container
////////////////////////////////////////////////////////////////////////////////
#endif
