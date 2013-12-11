#ifndef HOM3_MISC_STACK_VECTOR_HPP_
#define HOM3_MISC_STACK_VECTOR_HPP_
////////////////////////////////////////////////////////////////////////////////
/// \file \brief Stack allocated vector
////////////////////////////////////////////////////////////////////////////////
#include <vector>
#include "misc/stack_memory.hpp"
////////////////////////////////////////////////////////////////////////////////
namespace hom3 {
////////////////////////////////////////////////////////////////////////////////

namespace memory {

namespace stack {

/// \brief Stack allocated vector that can be returned from functions
/// \param T is the type of the elements stored in the vector
/// \param N is the max. number of elements that can be stored in the vector
/// \warning Because of the std::vector growth rules, \param I might need to be
/// larger than the actual maximum number of elements (this depends on the
/// implementation of std::vector)
template<class T, std::size_t N>
struct vector : std::vector<T, stack::allocator<T, N>> {
  using Base = std::vector<T, stack::allocator<T, N>>;
  using Other = std::vector<T>;

  /// \name Forward std::vector constructors
  ///@{
  inline explicit vector(typename Base::size_type count, const T& value)
    : Base(count, value, allocator())
  { (*this).reserve(N); }

  inline explicit vector() : vector(0, T{}) { (*this).reserve(N); }

  inline vector(const Other& other) : Base(other, allocator())
  { (*this).reserve(N); }
  inline vector(const vector& other) : vector(other) { (*this).reserve(N); }
  inline vector(vector&& other) : Base(std::move(other), allocator())
  { (*this).reserve(N); }
  inline vector(Other&& other) : Base(std::move(other), allocator())
  { (*this).reserve(N); }
  inline vector(std::initializer_list<T> init)
    : Base(std::move(init), allocator())
  { (*this).reserve(N); }
  ///@}

  using Base::operator=;  // Inherit std::vector assignment

private:
  stack::arena<T, N> memory_;
  inline stack::allocator<T, N> allocator()       { return {memory_}; }
  inline stack::allocator<T, N> allocator() const { return {memory_}; }
};

}  // namespace stack

}  // namespace memory

////////////////////////////////////////////////////////////////////////////////
}  // namespace hom3
////////////////////////////////////////////////////////////////////////////////
#endif
