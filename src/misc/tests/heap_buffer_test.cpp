/// \file Heap buffer tests
#include "gtest/gtest.h"
#include "globals.hpp"
////////////////////////////////////////////////////////////////////////////////
using namespace hom3;

/// \test Buffer usage example
TEST(heap_buffer_test, example_code) {
  // Create a memory buffer
  memory::Buffer buffer;
  // To use the buffer, create a buffer::Accesor with the required size:
  auto b = memory::buffer::Adapt<double>{buffer, 4};
  // Now you have a buffer \p that can hold 1000 doubles.
  // To access an element of the buffer you can use the () operator
  b(3) = 5.0;  // The () operator also performs out of bounds checking!
  // Buffer satisfy the container requirements, so you can iterate
  // over them using begin/end:
  for (auto&& i : b) { std::cout << i << "\n"; }
  // get a pointer to the raw buffer data:
  b.data();
  // as well as resizing the buffer, and asking the buffer size/capacity
  b.resize(10);
  ASSERT(b.size() == 10, "always true!");
  ASSERT(b.capacity() >= b.size(), "always true!");
}

template<class T>
void is_filled_with(memory::buffer::Adapt<T> b, T value) {
  Ind size = 0;
  for (auto&& i : b) {
    EXPECT_EQ(value, i);
    ++size;
  }
  EXPECT_EQ(b.size(), size);
}

void is_filled_with(memory::buffer::Adapt<Num> b, Num value) {
  Ind size = 0;
  for (auto&& i : b) {
    EXPECT_FLOAT_EQ(value, i);
    ++size;
  }
  EXPECT_EQ(b.size(), size);
}

template<class T> void basic_buffer_test(memory::Buffer& buf) {
  /// Create a buffer view of 10 Ts
  {
    memory::buffer::Adapt<T> b{buf, 10};
    EXPECT_EQ(b.capacity(), Ind{10});
    EXPECT_EQ(b.size(), Ind{10});
    b.fill();
    is_filled_with(b, T{});
    b.fill(T{3});
    is_filled_with(b, T{3});
  }

  /// The bufer is cleared up using RAII
  EXPECT_EQ(buf.capacity() / sizeof(T), std::size_t{10});
  EXPECT_EQ(buf.size(), std::size_t{0});

  /// A new view causes no new allocations
  {
    memory::buffer::Adapt<T> b{buf, 5};
    EXPECT_EQ(b.capacity(), Ind{10});  // capacity should still be 10 !
    EXPECT_EQ(b.size(), Ind{5});
  }
}

template<class T> void basic_buffer_test() {
  memory::Buffer b;
  basic_buffer_test<T>(b);
}

TEST(heap_buffer_test, basic_tests) {
  basic_buffer_test<Ind>();
  basic_buffer_test<SInd>();
  basic_buffer_test<Num>();
}

TEST(heap_buffer_test, buffer_reuse) {
  memory::Buffer buf;
  basic_buffer_test<Ind>(buf);
  basic_buffer_test<Num>(buf);
}
