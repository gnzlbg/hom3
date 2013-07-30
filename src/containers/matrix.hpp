#ifndef HOM3_CONTAINER_MATRIX_HPP_
#define HOM3_CONTAINER_MATRIX_HPP_
////////////////////////////////////////////////////////////////////////////////
/// Includes:
#include "../globals.hpp"
#include <type_traits>
////////////////////////////////////////////////////////////////////////////////

namespace container {

/// \todo rename tags/move them somewhere else/replace tag dispatch with crtp.

namespace matrix {

/// Column type traits
template<class T> struct traits {
  using value_type      = typename T::value_type;
  using reference       = typename T::reference;
  using const_reference = typename T::const_reference;
};

template<class T, int r, int c, int t, int p1, int p2>
struct traits<Eigen::Matrix<T,r,c,t,p1,p2>> {
 using value_type      = typename Eigen::Matrix<T,r,c,t,p1,p2>::Scalar;
 using reference       = typename Eigen::Matrix<T,r,c,t,p1,p2>::Scalar&;
 using const_reference = const reference;
};

namespace tag {
struct Cell {};
struct Node {};
struct NodeIndices {};
} // tag namespace

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
template<class C, class T, template <SInd> class V, SInd nd_ = 1> struct Matrix {

  /// #of columns
  inline static constexpr SInd nd(){ return nd_; }

  using This             = Matrix<C,T,V,nd()>;
  using data_container   = V<nd()>; ///< Underlying container type
  using reference        = typename matrix::traits<V<nd()>>::reference;
  using const_reference  = typename matrix::traits<V<nd()>>::const_reference;

  Matrix<C,T,V,nd()>() : c_(nullptr), data_(0,nd()), name_("unknown") {}
  Matrix<C,T,V,nd()>(const C* t, std::string name) : c_(t), data_(capacity_(),nd()), name_(name) {}
  Matrix<C,T,V,nd()>(const This& other) : c_(nullptr), data_(other()), name_(other.name()) {}
  Matrix<C,T,V,nd()>(This&& other) : Matrix<C,T,V,nd()>() { swap(*this, other); }
  This& operator=(This other) { swap(*this, other); return *this; }

  /// \brief Access element at index (\p row_i, \p col_j)
  inline reference       operator()(const Ind row_i,const SInd col_j = 0) {
    ASSERT(preconditions_(row_i,col_j),""); return data_(row_i,col_j);
  }
  inline const_reference operator()(const Ind row_i,const SInd col_j = 0) const {
    ASSERT(preconditions_(row_i,col_j),""); return const_reference(data_(row_i,col_j));
  }

  /// \brief Swaps data without swapping ownership:
  /// \todo swapping container ownership is easier and faster
  friend void swap(This& a, This& b) { std::swap(a.data_, b.data_); std::swap(a.name_,b.name_); }

  inline void init(const C* t) { c_ = t; data_.resize(capacity_(),nd()); }

  /// \brief Acces the underlying container type
  inline       data_container& operator()()       { return data_; }
  inline const data_container& operator()() const { return data_; }

  /// \brief Flips all bytes (work with containers of bools only!)
  inline void flip() { data_.flip(); } ///< flips all!
  inline void flip(const SInd d) { data_.flip(d); } //< flips all in dim nd_

  inline std::string name() const { return name_; }
 private:
  const C* c_;
  data_container data_;
  const std::string name_;

  inline std::string path() const { return name() + " in container " + c_->name(); }

  inline bool preconditions_(const Ind i, const SInd d) const {
    ASSERT(data_.size() != 0,"Uninitialized variables: " << path() << " has size = 0!");
    ASSERT(i < size_(),"Index " << i << " is out of bounds: " << path() << " has size = " << size_() << ".");
    ASSERT(i < capacity_(),"Index " << i << " in " << path() << " has capacity = " << capacity_() << ".");
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
};

// Alias (sugar) for array of node indices
template<class C> using NodeIndices = Matrix<C,matrix::tag::NodeIndices,IndM,1>;
////////////////////////////////////////////////////////////////////////////////
} // namespace cell
////////////////////////////////////////////////////////////////////////////////
#endif
