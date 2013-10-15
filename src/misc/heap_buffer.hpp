#ifndef HOM3_MISC_HEAP_BUFFER_HPP_
#define HOM3_MISC_HEAP_BUFFER_HPP_
////////////////////////////////////////////////////////////////////////////////
/// \file \brief A heap allocated buffer
///
/// For an example, read TEST(heap_buffer_test, example_code) in
/// misc/tests/heap_buffer_test.cpp.
///
////////////////////////////////////////////////////////////////////////////////
/// Include:
#include <vector>
#include "globals.hpp"
////////////////////////////////////////////////////////////////////////////////
namespace hom3 { namespace memory {
////////////////////////////////////////////////////////////////////////////////

using Buffer = std::vector<char>;

namespace buffer {

/// \brief Adapt's a memory buffer to the type T
template <class T> struct Adapt {
  using This            = Adapt<T>;
  using value_type      = T;
  using iterator        = T*;
  using const_iterator  = T const *;
  using reference       = T&;
  using const_reference = T const&;

  Adapt(Buffer& buffer, const Ind noTs = 0) : buf_(buffer) { resize(noTs); }
  ~Adapt() { buf_.clear(); }
  Adapt()                = delete;

  /// \brief Buffer capacity to hold T's
  inline Ind capacity() const noexcept { return buf_.capacity() / sizeof(T); }
  /// \brief Makes sure buffer is able to hold \p noTs \p T 's.
  inline void resize(const Ind noTs) { buf_.resize(noTs * sizeof(T)); }
  /// \brief Makes sure buffer is able to hold \p noTs \p T 's.
  inline Ind size() const noexcept { return buf_.size() / sizeof(T); }
  /// \brief Access \p i 'th T element of the buffer
  inline reference operator()(const Ind i) noexcept {
    ASSERT(i >= 0 && i < size(), "Element i = " << i << " is out of bounds!");
    return reinterpret_cast<reference>(buf_[i * sizeof(T)]);
  }
  inline value_type operator()(const Ind i) const noexcept {
    ASSERT(i >= 0 && i < size(), "Element i = " << i << " is out of bounds!");
    return reinterpret_cast<reference>(buf_[i * sizeof(T)]);
  }

  /// \brief Fill
  inline void fill(T value = T{}) noexcept { boost::fill(*this, value); }

  /// \brief Pointer to raw data
  inline iterator data() noexcept
  { return reinterpret_cast<iterator>(buf_.data()); }
  inline const_iterator data() const noexcept
  { return reinterpret_cast<iterator>(buf_.data()); }

  inline iterator begin() noexcept { return data(); }
  inline iterator end() noexcept { return data() + size(); }
  inline const_iterator begin() const noexcept { return data(); }
  inline const_iterator end() const noexcept { return data() + size(); }

 private:
  Buffer& buf_;
};

}  // namespace buffer

////////////////////////////////////////////////////////////////////////////////
}  // namespace memory
}  // namespace hom3
////////////////////////////////////////////////////////////////////////////////
#endif
