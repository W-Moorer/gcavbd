#include <gtest/gtest.h>

#include "gcavbd/core/types.h"

TEST(MathTest, OnePlusOne) {
  EXPECT_EQ(1 + 1, 2);
}

TEST(MathTest, Vec2Addition) {
  gcavbd::Vec2 a(1.0, 2.0);
  gcavbd::Vec2 b(3.0, 4.0);
  gcavbd::Vec2 c = a + b;

  EXPECT_DOUBLE_EQ(c(0), 4.0);
  EXPECT_DOUBLE_EQ(c(1), 6.0);
}
