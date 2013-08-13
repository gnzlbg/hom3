#ifndef HOM3_MEMORY_STACK_HPP_
#define HOM3_MEMORY_STACK_HPP_
////////////////////////////////////////////////////////////////////////////////
/// \file \brief Classes and helpers to manage memory on the stack
///
/// Example: create std::list<int> and std::vector<int> using stack memory
/// to allocate a maximum of 20 elements each:
/// \code
/// memory::stack::arena<int, 20> listMemory;
/// auto smallList = memory::stack::make<std::list>(listMemory);
/// memory::stack::arena<int, 20> vectorMemory;
/// auto smallVector = memory::stack::make<std::vector>(vectorMemory);
/// \endcode
///
/// Copyright: Howard Hinnant
/// - http://home.roadrunner.com/~hinnant/short_alloc.h
/// - http://home.roadrunner.com/~hinnant/stack_alloc.html
///
/// \warning Modifications to Howard's code
/// - <cassert> has been replaced by HOM3::ASSERT
/// - code has been reformated
/// - short_alloc renamed to allocator
/// - added make function
/// - added type T to arena, recompute N as N = sizeof(T) * M
/// - terminate if runs out of stack memory
////////////////////////////////////////////////////////////////////////////////
namespace hom3 {
////////////////////////////////////////////////////////////////////////////////
namespace memory {
////////////////////////////////////////////////////////////////////////////////
namespace stack {
////////////////////////////////////////////////////////////////////////////////

template <class T, std::size_t M>
struct arena {
  static const std::size_t N = sizeof(T) * M;
  arena() noexcept : ptr_(buf_) {}
  ~arena() {ptr_ = nullptr;}
  arena(const arena&) = delete;
  arena& operator=(const arena&) = delete;

  char* allocate(std::size_t n);
  void deallocate(char* p, std::size_t n) noexcept;

  static constexpr std::size_t size() { return N; }
  std::size_t used() const { return static_cast<std::size_t>(ptr_ - buf_); }
  void reset() { ptr_ = buf_; }

 private:
  static const std::size_t alignment = 16;
  alignas(alignment) char buf_[N];
  char* ptr_;

  std::size_t align_up(std::size_t n) noexcept
  { return n + (alignment-1) & ~(alignment-1); }

  bool pointer_in_buffer(char* p) noexcept
  { return buf_ <= p && p <= buf_ + N; }
};

template <class T, std::size_t M>
char* arena<T,M>::allocate(std::size_t n) {
  ASSERT(pointer_in_buffer(ptr_), "allocator has outlived arena");
  n = align_up(n);
  if (buf_ + N - ptr_ >= static_cast<std::ptrdiff_t>(n)) {
    char *const r = ptr_;
    ptr_ += n;
    return r;
  }
  TERMINATE("not enough memory " + AT_);
  //return static_cast<char*>(::operator new(n));
}

template <class T, std::size_t M>
void arena<T,M>::deallocate(char* p, std::size_t n) noexcept {
  ASSERT(pointer_in_buffer(ptr_), "allocator has outlived arena");
  if (pointer_in_buffer(p)) {
    n = align_up(n);
    if (p + n == ptr_) {
      ptr_ = p;
    }
  } else {
    TERMINATE("should never happen");
    //::operator delete(p);
  }
}

template <class T, std::size_t N>
struct allocator {
  typedef T value_type;

  template <class _Up>
  struct rebind {typedef allocator<_Up, N> other;};

  allocator(arena<T,N>& a) noexcept : a_(a) {}
  template <class U>
  allocator(const allocator<U, N>& a) noexcept : a_(a.a_) {}
  allocator(const allocator&) = default;
  allocator& operator=(const allocator&) = delete;

  T* allocate(std::size_t n)
  { return reinterpret_cast<T*>(a_.allocate(n*sizeof(T))); }
  void deallocate(T* p, std::size_t n) noexcept
  { a_.deallocate(reinterpret_cast<char*>(p), n*sizeof(T)); }

  template <class T1, std::size_t N1, class U, std::size_t M>
  friend
  bool
  operator==(const allocator<T1, N1>& x, const allocator<U, M>& y) noexcept;

  template <class U, std::size_t M> friend struct allocator;

 private:
  arena<T,N>& a_;
};

template <class T, std::size_t N, class U, std::size_t M>
inline
bool operator==(const allocator<T, N>& x, const allocator<U, M>& y) noexcept
{ return N == M && &x.a_ == &y.a_; }

template <class T, std::size_t N, class U, std::size_t M>
inline
bool operator!=(const allocator<T, N>& x, const allocator<U, M>& y) noexcept
{ return !(x == y); }

template<template<class,class> class C, class T, std::size_t N>
using Container = C<T,allocator<T,N>>;

template<template<class,class> class C, class T, std::size_t N>
Container<C,T,N> make(arena<T,N>& arena_) {
  return Container<C,T,N>{allocator<T,N>{arena_}};
}


////////////////////////////////////////////////////////////////////////////////
} // namespace stack
////////////////////////////////////////////////////////////////////////////////
} // namespace memory
////////////////////////////////////////////////////////////////////////////////
} // namespace hom3
////////////////////////////////////////////////////////////////////////////////
#endif
