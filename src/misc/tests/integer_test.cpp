/// \file Tests the Integer wrapper
#include "gtest/gtest.h"
#include "globals.hpp"
////////////////////////////////////////////////////////////////////////////////
using namespace hom3;

TEST(types_test, integer_conversions) {
  /// Should be constructible from its underlying type
  auto a = Integer<Ind>{2};

  /// Should be explicitly convertible to its underlying type
  Ind b = a();
  Ind c = 2;
  EXPECT_EQ(b, c);
  auto lambda = [](Ind i) { return i; };
  EXPECT_EQ(lambda(a()), c);

  /// Should not be implicitly convertible to its underlying type
  // Ind d = a; // does not compile :)
  // lambda(a); // does not compile :)

  /// Should not implicitly convert to integer types of different families
  struct class1 {}; struct class2 {};
  auto e = Integer<Ind, class1>{2};
  auto f = Integer<Ind, class2>{3};
  // e = f; // does not compile :)
  e = f();
}

TEST(types_test, integer_primitive_cast) {
  struct class1 {};
  using T1 = Integer<Ind, class1>;
  struct class2 {};
  using T2 = Integer<Ind, class2>;

  auto a = T1{2};
  auto b = T2{3};

  primitive_cast(a) = primitive_cast(b);
  EXPECT_EQ(a, T1{3});
}

template<class T> void test_integer_compound_assignment() {
  auto i1 = T{1};
  auto i2 = T{2};
  i1 += i1;
  EXPECT_EQ(i1, i2);
  i2 -= i2;
  EXPECT_EQ(i2, T{0});
  i1 *= T{2};
  EXPECT_EQ(i1, T{4});
  i1 /= T{2};
  EXPECT_EQ(i1, T{2});
}

TEST(types_test, integer_compound_assignment) {
  test_integer_compound_assignment<Integer<Ind>>();
  test_integer_compound_assignment<Integer<Int>>();
}

template<class T> void test_unsigned_integer_arithmetic() {
  auto i1 = T{1};
  auto i2 = i1 + i1 + i1;
  auto i3 = i1 * T{3};
  EXPECT_EQ(i2, i3);
  EXPECT_EQ(i3 / i2, T{1});
  auto i4 = i3 - T{2} * i1;
  EXPECT_EQ(i4, i1);
}

TEST(types_test, integer_arithmetic) {
  test_unsigned_integer_arithmetic<Integer<Ind>>();
  test_unsigned_integer_arithmetic<Integer<Int>>();

  // signed integers:
  using I = Integer<Int>;
  auto i1 = I{1};
  auto i3 = i1 - I{2} * i1;
  EXPECT_EQ(i3, -i1);
}

template<class T> void test_unsigned_integer_increment_operators() {
  auto i1 = T{1};
  auto i2 = T{2};
  ++i1;
  EXPECT_EQ(i1, i2);
  auto i3 = i2;
  EXPECT_EQ(i3, i2++);
  EXPECT_EQ(i2, T{3});
  --i2;
  EXPECT_EQ(i2, i1);
  auto i4 = i1;
  EXPECT_EQ(i4, i1--);
}

TEST(types_test, integer_increment_operators) {
  test_unsigned_integer_increment_operators<Integer<Ind>>();
  test_unsigned_integer_increment_operators<Integer<Int>>();
}

template<class T> void test_comparison_operators() {
  auto a = T{2};
  auto b = T{3};
  auto c = T{5};
  auto d = T{1};
  EXPECT_TRUE(a + b == c);
  EXPECT_TRUE(a < b);
  EXPECT_TRUE(a <= b);
  EXPECT_TRUE(b <= b);
  EXPECT_TRUE(b >= b);
  EXPECT_TRUE(b > a);
  EXPECT_TRUE(b >= a);
  EXPECT_TRUE(b == b);
  EXPECT_TRUE(b != a);
  EXPECT_TRUE(b - a == d);
  EXPECT_TRUE(a * b == c + d);

  auto a_old = a;
  a += b;
  EXPECT_TRUE(a == c);
  b -= a_old;
  EXPECT_TRUE(b == d);
}

TEST(types_test, integer_comparison_operators) {
  test_comparison_operators<Integer<Ind>>();
  test_comparison_operators<Integer<Int>>();
}
