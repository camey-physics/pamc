#include <gtest/gtest.h>
#include <memory>
#include <cstddef>

#include "MemoryPool.hpp"

TEST(MemoryPoolTest, AllocateSingleDoubles) {
  const std::size_t pool_size = 10;
  MemoryPool<double> pool(pool_size);

  double* ptrs[pool_size];

  // Allocate 10 separate doubles
  for (std::size_t i = 0; i < pool_size; ++i) {
    ptrs[i] = pool.allocate(1);
    EXPECT_NE(ptrs[i], nullptr);
    if (i > 0) {
      // Should point to the next double
      EXPECT_EQ(ptrs[i], ptrs[i - 1] + 1);
    }
  }

  // Pool should now be full
  EXPECT_EQ(pool.size(), pool_size);
  EXPECT_EQ(pool.capacity(), pool_size);
}

TEST(MemoryPoolTest, AllocateEightDoubles) {
  const std::size_t pool_size = 10;
  MemoryPool<double> pool(pool_size *8);

  double* ptrs[pool_size];

  // Allocate 10 separate doubles
  for (std::size_t i = 0; i < pool_size; ++i) {
    ptrs[i] = pool.allocate(8);
    EXPECT_NE(ptrs[i], nullptr);
    if (i > 0) {
      // Should point to the next double
      EXPECT_EQ(ptrs[i], ptrs[i - 1] + 8);
    }
  }

  // Pool should now be full
  EXPECT_EQ(pool.size(), pool_size *8);
  EXPECT_EQ(pool.capacity(), pool_size *8);
}

TEST(MemoryPoolTest, ResetAllowsReuse) {
  MemoryPool<int> pool(10);
  int* a = pool.allocate(10);
  pool.reset();
  int* b = pool.allocate(10);
  // Memory should be reused
  EXPECT_EQ(a, b);
}

TEST(MemoryPoolTest, OverAllocateTriggersAssert) {
  MemoryPool<int> pool(5);
  pool.allocate(5);
  EXPECT_DEATH(pool.allocate(1), ".*");
}