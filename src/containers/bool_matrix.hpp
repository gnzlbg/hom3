#ifndef HOM3_CONTAINERS_BOOL_MATRIX_HPP_
#define HOM3_CONTAINERS_BOOL_MATRIX_HPP_
////////////////////////////////////////////////////////////////////////////////
/// \file \brief Bool matrix class
////////////////////////////////////////////////////////////////////////////////

namespace container {

/// \brief Bool Matrix type: ColMajor
template <SInd nd> struct BoolMatrix {

  using This            = BoolMatrix<nd>;
  using container       = std::vector<bool>;
  using value_type      = typename container::value_type;
  using reference       = typename container::reference;
  using const_reference = typename container::const_reference;

  BoolMatrix<nd>(const Ind n, const SInd /*nd_*/) {
    //ASSERT(nd_ == nd, "Wrong number of dimensions!");
    for(auto&& v : data_) {
      v.resize(n);
    }
  };

  BoolMatrix<nd>(const This& other) : data_(other.data_) {}
  This& operator=(This other) { swap(*this, other); return *this; }
  friend void swap(This& a, This& b) { std::swap(a.data_,b.data_); }

  inline reference       operator()(const Ind i,const SInd d = 0)       {
    return data_[d][i];
  }
  inline const_reference operator()(const Ind i,const SInd d = 0) const {
    return data_[d][i];
  }

  inline bool check_size() const {
    ASSERT(data_.size() > 0, "Can't check size on empty data_");
    const auto size_ = data_[0].size();
    for(auto&& v : data_) {
      if( v.size() != size_ ) {
        return false;
      }
    }
    return true;
  }

  std::size_t size() const {
    ASSERT( check_size(), "Inconsistent size");
    return nd * data_[0].size();
  }

  void resize(const Ind n, const SInd nd_) {
    ASSERT(nd_ == nd, "Wrong number of dimensions!");
    for(auto&& v : data_) {
      v.resize(n);
    }
  }

 private:
  std::array<container,nd> data_;
};

} // container namespace

using container::BoolMatrix;

////////////////////////////////////////////////////////////////////////////////
#endif
