#include "models/IsingModel.hpp"
#include "models/Ising3DHelpers.hpp"
#include "MemoryBlock.hpp"
#include "MemoryPool.hpp"
#include <gtest/gtest.h>
#include <typeindex>

class TestIsingModelPool : public ::testing::Test {
 protected:
  int L = 5;
  int num_spins = L * L * L;
  int num_neighbors = 6;
  double J = 1.0;
  std::vector<int> neighbor_table = initializeNeighborTable3D(L);
  std::vector<double> bond_table = std::vector<double>(num_spins * num_neighbors, J);
  SharedModelData<IsingModel> shared_data = SharedModelData<IsingModel>(
      L, num_spins, num_neighbors, neighbor_table.data(), bond_table.data());
};

TEST_F(TestIsingModelPool, ReportsStorageRequirements) {
  IsingModel model(shared_data);
  std::size_t expected_spins = L * L * L;

  auto reqs = model.storageRequirements(L);
  ASSERT_EQ(reqs.size(), 1);
  EXPECT_EQ(reqs[0].type, std::type_index(typeid(int)));
  EXPECT_EQ(reqs[0].count, expected_spins);
}

TEST_F(TestIsingModelPool, Construct) {
  MemoryPool<int> pool = MemoryPool<int>(num_spins);
  IsingModel model(shared_data, pool.allocate(num_spins));

  ASSERT_TRUE(model.usesExternalPool());

  for (int i = 0; i < num_spins; ++i) {
    EXPECT_EQ(model.getSpin(i), 1);
  }
}

TEST_F(TestIsingModelPool, CopyState) {
  MemoryPool<int> pool = MemoryPool<int>(num_spins *2);
  IsingModel model(shared_data, pool.allocate(num_spins));
  IsingModel model2(shared_data, pool.allocate(num_spins));
  gsl_rng* r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r, 42);
  model.initializeState(r);
  model2.copyStateFrom(model);

  for (int i = 0; i < num_spins; ++i) {
    EXPECT_EQ(model.getSpin(i), model2.getSpin(i));
  }

  EXPECT_NE(model.getSpins_(), model2.getSpins_());

  gsl_rng_free(r);
}
