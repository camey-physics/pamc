#include <gtest/gtest.h>

#include "SharedModelData.hpp"
#include "models/IsingModel.hpp"

IsingModel createDummyModel() {
  int dummy_neighbors[6] = {0};
  double dummy_bonds[6] = {1.0};
  SharedModelData<IsingModel> d(1, 1, 6, dummy_neighbors, dummy_bonds);
  return IsingModel(d);
}

TEST(IsingModelTest, Construct) {
  const int L = 4;
  const int num_spins = L * L * L;
  const int num_neighbors = 6;

  // Create dummy neighbor and bond tables, fill with dummy data
  int neighbor_table[num_spins * num_neighbors] = {0};
  double bond_table[num_spins * num_neighbors] = {0.0};

  SharedModelData<IsingModel> data(L, num_spins, num_neighbors, neighbor_table,
                                   bond_table);
  IsingModel model(data);

  // Check default initialization of spins_
  for (int i = 0; i < num_spins; ++i) {
    EXPECT_EQ(model.getSpin(i), 1);
  }
}
