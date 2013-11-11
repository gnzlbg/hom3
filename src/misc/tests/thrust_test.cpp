/// \file \brief Implements tests for Thrust related functionality
#include "misc/test.hpp"
#include "globals.hpp"
#include "misc/compute.hpp"
#include "misc/meta.hpp"

#include <vector>
////////////////////////////////////////////////////////////////////////////////
using namespace hom3;

int main(int argc, char* argv[]) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

/// Note: zipped is an implementation detail of iterator: skip it on a first
/// read of the class, the public interface comes right afterwards.
template<class column_tags, class column_types> struct zipped {
  /// \brief (column_tags, column_types) -> (column_tags, P(column_types))
  template<template <class> class P>
  using transformed = typename fusion::result_of::as_map<
    typename mpl::transform<
      typename mpl::transform<column_types, P<mpl::_1>>::type,
      column_tags, fusion::pair<mpl::_2, mpl::_1>
        >::type
    >::type;

    /// \brief (column_tags, column_types) -> (column_tags, P(column_types))
  template<template <class> class P>
  using transformed = typename fusion::result_of::as_map<
    typename mpl::transform<
      typename mpl::transform<column_types, P<mpl::_1>>::type,
      column_tags, fusion::pair<mpl::_2, mpl::_1>
        >::type
    >::type;

  using value_map = typename fusion::result_of::as_map<
    typename mpl::transform<
      column_types, column_tags, fusion::pair<mpl::_2, mpl::_1>
      >::type
    >::type;

  using reference_map = transformed<boost::add_reference>;
};

template<class column_types, class column_tags,
         template <class> class iterator_type>
struct iterator {
  using col_types = column_types;
  using col_tags  = column_tags;

  template<class T> struct iterator_lambda {  using type = iterator_type<T>; };
  using iterator_map = typename zipped<col_tags, col_types>::template
                         transformed<iterator_lambda>;

  /// \name Iterator traits
  ///@{
  using value_type    = typename zipped<col_tags, col_types>::value_map;
  using reference_map = typename zipped<col_tags, col_types>::reference_map;
  using iterator_category = std::random_access_iterator_tag;
  using difference_type   = Int;
  using pointer = iterator_map;
  ///@}

  using This = iterator<col_types, col_tags, iterator_type>;

  struct reference : reference_map {
    reference(reference& r) noexcept : reference_map(r) {};
    reference(value_type& v) noexcept : reference_map(ref_from_val(v)) {}
    template<class... Ts> reference(Ts&&... ts) noexcept
      : reference_map(std::forward<Ts>(ts)...) {}

    template<class T> reference& operator=(T&& other) noexcept {
      fusion::for_each(other, fusion::map_assigner<reference>(*this));
      return *this;
    }
  };

  iterator() = default;
  template<class... Ts>
  explicit iterator(Ts... ts) noexcept : it_(std::move(ts)...) {}
  explicit iterator(iterator_map it) noexcept : it_(std::move(it)) {}

  /// \name Comparison operators (==, !=, <, >, <=, >=)
  ///@{
  friend inline bool operator==(const This& l, const This& r) noexcept
  { return fusion::all(fusion::zip(l.it_, r.it_), Eq()); }
  friend inline bool operator<=(const This& l, const This& r) noexcept
  { return fusion::all(fusion::zip(l.it_, r.it_), LEq()); }
  friend inline bool operator>=(const This& l, const This& r) noexcept
  { return fusion::all(fusion::zip(l.it_, r.it_), GEq()); }

  friend inline bool operator!=(const This& l, const This& r) noexcept
  { return !(l == r); }
  friend inline bool operator< (const This& l, const This& r) noexcept
  { return !(l >= r); }
  friend inline bool operator> (const This& l, const This& r) noexcept
  { return !(l <= r); }
  ///@}

  /// \name Traversal operators (++, +=, +, --, -=, -)
  ///@{
  inline This& operator++() noexcept   {
    it_ = fusion::transform(it_, AdvFwd());
    return *this;
  }
  inline This operator++(int) noexcept {
    const auto it = *this;
    ++(*this);
    return it;
  }
  inline This& operator+=(const difference_type value) noexcept {
    it_ = fusion::transform(it_, Increment(value));
    return *this;
  }
  inline This operator--(int) noexcept {
    const auto it = *this;
    --(*this);
    return it;
  }
  inline This& operator-=(const difference_type value) noexcept {
    it_ = fusion::transform(it_, Decrement(value));
    return *this;
  }
  inline This& operator--()   noexcept {
    it_ = fusion::transform(it_, AdvBwd());
    return *this;
  }

  friend inline This operator-(This a, const difference_type value) noexcept
  { return a -= value; }
  friend inline This operator+(This a, const difference_type value) noexcept
  { return a += value; }

  friend inline
  difference_type operator-(const This l, const This r) noexcept {
    const auto result
      = fusion::at_c<0>(l.it_).second - fusion::at_c<0>(r.it_).second;
    ASSERT(fusion::accumulate(fusion::zip(l.it_, r.it_), result,
                              check_distance()),
           "column types have different distances!");
    return result;
  }
  ///@}

  /// \name Access operators
  ///@{
  inline reference operator* ()       noexcept { return ref_from_it(it_); }
  inline reference operator* () const noexcept { return ref_from_it(it_); }
  inline reference operator->()       noexcept { return *(*this);         }
  inline reference operator[](const difference_type v) noexcept
  { return *fusion::transform(it_, Increment(v)); }
  ///@}

  /// \brief Makes an iterator from an iterator sequence.
  ///
  /// \warning The iterator sequence must be in the same order and have the same
  /// types as the iterator column_types (this is checked at compile time).
  template<class... Ts> static inline This make(Ts&&... ts) noexcept {
    using types = mpl::vector<std::remove_reference_t<Ts>...>;
    using iterator_types = typename mpl::transform<
      column_types, iterator_lambda<mpl::_1>>::type;
    static_assert(mpl::equal<types, iterator_types>::value,
                  "the types of the iterators provided do not match with the \
                  column type col_type of this iterator");
    return This{make_iterator_map(std::forward<Ts>(ts)...)};
  }

 private:
  iterator_map it_;  ///< Tuple of iterators over column types

  /// \name Make iterator implementation

  template<class... Ts>
  static inline iterator_map make_iterator_map(Ts&&... ts) noexcept
  { return fusion::as_map(make_iterator_tuple(std::forward<Ts>(ts)...)); }

  struct zip2pair {
    template<class T> inline auto operator()(T&& i) const noexcept {
      using key = std::remove_reference_t<decltype(fusion::at_c<0>(i))>;
      auto value = fusion::at_c<1>(std::forward<T>(i));
      return fusion::make_pair<key>(value);
    }
  };

  template<class... Ts>
  static inline auto make_iterator_tuple(Ts&&... ts) noexcept {
    const auto fusion_vector_val = fusion::make_vector(std::forward<Ts>(ts)...);
    return fusion::as_vector
        (fusion::transform
         (fusion::zip(col_tags(), fusion_vector_val), zip2pair()));
  }
  ///@}

  /// \name Wrappers
  ///@{

  /// \brief Equality
  struct Eq { template<class T>   inline bool operator()(T&& i) const noexcept
  { return fusion::at_c<0>(i).second == fusion::at_c<1>(i).second; }};

  /// \brief Less Equal Than
  struct LEq { template<class T>  inline bool operator()(T&& i) const noexcept
  { return fusion::at_c<0>(i).second <= fusion::at_c<1>(i).second; }};

  /// \brief Greater Equal Than
  struct GEq { template<class T>  inline bool operator()(T&& i) const noexcept
  { return fusion::at_c<0>(i).second >= fusion::at_c<1>(i).second; }};

  /// \brief Advance forward
  struct AdvFwd { template<class T> inline auto operator()(T i) const noexcept
  { return i.second = ++(i.second); }};

  /// \brief Advance backward
  struct AdvBwd { template<class T> inline auto operator()(T i) const noexcept
  { return i.second = --(i.second); }};

  /// \brief Increment by value
  struct Increment {
    Increment(const difference_type value) noexcept : value_(value) {}
    template<class T>
    inline std::remove_reference_t<T> operator()(T&& i) const noexcept
    { return {i.second + value_}; }
    const difference_type value_;
  };

  /// \brief Decrent by value
  struct Decrement {
    Decrement(const difference_type value) : value_(value) {}
    template<class T>
    inline std::remove_reference_t<T> operator()(T&& i) const noexcept
    { return {i.second - value_}; }
    const difference_type value_;
  };

  struct check_distance {
    template<class T, class U>
    inline T operator()(const T acc, const U i) const noexcept {
      const auto tmp = fusion::at_c<0>(i).second
                       - fusion::at_c<1>(i).second;
      ASSERT(acc == tmp, "column types have different sizes!!");
      return tmp;
    }
  };

  ///@}

  template<class U> struct val_to_ref {
    val_to_ref(U& v_) : v(v_) {}
    template<class T>
    inline auto operator()(const T&) const noexcept {
      using key = typename T::first_type;
      using value = typename T::second_type;
      return fusion::make_pair<key, std::add_lvalue_reference_t<value>>
          (fusion::at_key<key>(v));
    }
    U& v;
  };

  template<class U> struct it_to_ref {
    it_to_ref(U v_) : i(v_) {}
    template<class T>
    inline auto operator()(const T&) const noexcept {
      using key = typename T::first_type;
      using value
        = std::remove_reference_t<decltype(*typename T::second_type())>;
      return fusion::make_pair<key, std::add_lvalue_reference_t<value>>
          (*fusion::at_key<key>(i));
    }
    U i;
  };

  static inline reference_map ref_from_val(value_type& v)          noexcept
  { return fusion::transform(v, val_to_ref<value_type>{v}); }
  static inline reference_map ref_from_it(iterator_map& v)       noexcept
  { return fusion::transform(v, it_to_ref<iterator_map>{v}); }
  static inline reference_map ref_from_it(const iterator_map& v) noexcept
  { return fusion::transform(v, it_to_ref<iterator_map>{v}); }
};


struct X {}; struct Y {}; struct I {}; struct B {};

template<class T>
struct HDVector {
  host::vector<T> host_;
  device::vector<T> device_;
};

template<class T> HDVector<T> make_hdvector(int size, T v = T{}) {
  HDVector<T> tmp;
  tmp.host_.resize(size, v);
  tmp.device_.resize(size, v);
  return tmp;
}

namespace test {

struct Object {
  float x;
  float y;
  int   i;
  bool  b;
};

}

namespace keys {

struct x {}; struct y{}; struct i{}; struct b{};

}

BOOST_FUSION_ADAPT_ASSOC_STRUCT(
    test::Object,
    (float, x, keys::x)
    (float, y, keys::y)
    (int,   i, keys::i)
    (bool,  b, keys::b)
    )

template<template <class> class UnderlyingContainer, class Object> struct StructOfArrays {
  template<class O> struct to_key { using type = typename O::first_type; };
  template<class O> struct to_val { using type = typename O::second_type; };

  template<class V> struct container_of { using type = UnderlyingContainer<V>; };

  template<class O> struct container_it
  { using type = typename container_of<O>::iterator; };

  using arrays = typename fusion::result_of::as_map<
    typename mpl::transform<
      typename mpl::transform<Object, typename container_of<typename to_val<mpl::_1>::type>::type,
                              typename mpl::transform<Object, to_key<mpl::_1>>::type,
                              fusion::pair<mpl::_2, mpl::_1>
                              >::type
      >::type
    >::type;

  using keys   = typename mpl::transform<Object, to_key<mpl::_1>>::type;
  using values = typename mpl::transform<Object, to_val<mpl::_1>>::type;

  using it = iterator<values, keys, container_it>;

  inline it begin() {
    return it::make(fusion::transform(data_, [](auto i){ return i.begin(); }));
  }
  inline it end()   {
    return it::make(fusion::transform(data_, [](auto i){ return i.end(); }));
  }

 private:
  arrays data_;
};

struct Container {
  Container()
    : X_(make_hdvector<float>(10, 1.))
    , Y_(make_hdvector<float>(10))
    , I_(make_hdvector<int>  (10))
    , B_(make_hdvector<bool> (10))
  {}

  HDVector<float> X_;
  HDVector<float> Y_;
  HDVector<int>   I_;
  HDVector<bool>  B_;

  using column_types = mpl::vector<float, float, int, bool>;
  using column_tags  = mpl::vector<    X,     Y,   I,    B>;

  template<class T> using hv_it = typename host::vector<T>::iterator;
  template<class T> using dv_it = typename device::vector<T>::iterator;

  using host_iterator   = iterator<column_types, column_tags, hv_it>;
  using device_iterator = iterator<column_types, column_tags, dv_it>;

  inline host_iterator host_begin() noexcept {
    return host_iterator::make(X_.host_.begin(), Y_.host_.begin(),
                               I_.host_.begin(), B_.host_.begin());
  }

  inline host_iterator host_end() noexcept {
    return host_iterator::make(X_.host_.end(), Y_.host_.end(),
                               I_.host_.end(), B_.host_.end());
  }

  inline device_iterator device_begin() noexcept {
    return device_iterator::make(X_.device_.begin(), Y_.device_.begin(),
                                 I_.device_.begin(), B_.device_.begin());
  }

  inline device_iterator device_end() noexcept {
    return device_iterator::make(X_.device_.end(), Y_.device_.end(),
                                 I_.device_.end(), B_.device_.end());
  }

  inline auto host() noexcept
  { return boost::make_iterator_range(host_begin(), host_end()); }

  inline auto device() noexcept
  { return boost::make_iterator_range(device_begin(), device_end()); }
};


struct sort_cmp {
  template<class T, class U>
  bool operator()(T&& t, U&& u) const {
    return fusion::at_key<X>(t) > fusion::at_key<X>(u);
  }
};

TEST(thrust_test, initialization) {
  Container c;
  auto it = c.host_begin();
  auto first_elements = *it;

  typename Container::host_iterator::value_type value;
  fusion::at_key<X>(value) = 3.;
  fusion::at_key<Y>(value) = 3.;
  fusion::at_key<I>(value) = 3;
  fusion::at_key<B>(value) = true;

  *it = value;

  int count = 0;
  for (auto i : c.host()) {
    fusion::at_key<X>(i) = double(count++);
    std::cout << "("  << fusion::at_key<X>(i)
              << ", " << fusion::at_key<Y>(i)
              << ", " << fusion::at_key<I>(i)
              << ", " << fusion::at_key<B>(i) << ")\n";
  }

  algorithm::stable_sort(c.host(), sort_cmp());

  std::cout << "SORTED:\n";
  for (auto i : c.host()) {
    std::cout << "("  << fusion::at_key<X>(i)
              << ", " << fusion::at_key<Y>(i)
              << ", " << fusion::at_key<I>(i)
              << ", " << fusion::at_key<B>(i) << ")\n";
  }

  device::vector<float> X(9), Y(9);

  compute::sequence(X.begin(), X.end());
  compute::sequence(Y.begin(), Y.end());
  const auto a = 1.0;

  std::cout << "Vectors before:\n";
  for(auto i : zip(X, Y)) {
    std::cout << "(" << boost::get<0>(i) << ", " << boost::get<1>(i) << ")\n";
  }



  /// SAXPY
  thrust::transform(X.begin(), X.end(), Y.begin(), Y.begin(),
                    [=](const float& x, const float& y) { return a * x + y; });

  std::cout << "Vectors after:\n";
  for(auto i : zip(X, Y)) {
    std::cout << "(" << boost::get<0>(i) << ", " << boost::get<1>(i) << ")\n";
  }
}
