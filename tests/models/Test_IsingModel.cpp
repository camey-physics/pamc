#include <gtest/gtest.h>
#include <gsl/gsl_rng.h>

#include <vector>
#include "SharedModelData.hpp"
#include "models/IsingModel.hpp"
#include "models/Ising3DHelpers.hpp"

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

TEST(IsingModelTest, CreateNeighborTable) {
  int L = 4;
  std::vector<int> table = initializeNeighborTable3D(L);
  EXPECT_EQ(table.size(), L * L * L * 6);
}

// Test the neighbors of spin 0 and L *L
TEST(IsingModelTest, NeighborTable) {
  int L = 4;
  int spin = 0;
  int num_neighbors = 6;
  std::vector<int> table = initializeNeighborTable3D(L);
  std::vector<int> expectedNeighbors = {(L-1)*L*L, L*L, (L-1)*L, L, L-1, 1};
  for (int n = 0; n < num_neighbors; ++n) {
    EXPECT_EQ(table[spin *num_neighbors + n], expectedNeighbors[n]);
  }

  expectedNeighbors = {0, 32, 28, 20, 19, 17};
  spin = L *L;
    for (int n = 0; n < num_neighbors; ++n) {
    EXPECT_EQ(table[spin *num_neighbors + n], expectedNeighbors[n]);
  }

}

TEST(IsingModelTest, CalcEnergy) {
  const int L = 5;
  const int num_spins = L * L * L;
  const int num_neighbors = 6;

  std::vector<int> neighbor_table = initializeNeighborTable3D(L);
  std::vector<double> bond_table(num_spins * num_neighbors, -1.0);

  SharedModelData<IsingModel> data(L, num_spins, num_neighbors, neighbor_table.data(),
                                   bond_table.data());

  // Calculate energy of all spins up configuration
  IsingModel model(data);
  EXPECT_NEAR(model.calcEnergy(), -3.0, 1e-10);
  // Calculate energy with one spin down
  int ind = index3D(1, 0, 0, L);
  model.setSpin(ind, -1);
  EXPECT_NEAR(model.calcEnergy(), -3.0 + 12.0 /num_spins, 1e-10);
  // Calculate energy with neighboring spin down
  ind = index3D(0, 0, 0, L);
  model.setSpin(ind, -1);
  EXPECT_NEAR(model.calcEnergy(), -3.0 + 20.0 /num_spins, 1e-10);
  // Calculate energy with second neighboring spin down
  ind = index3D(0, L - 1, 0, L);
  model.setSpin(ind, -1);
  EXPECT_NEAR(model.calcEnergy(), -3.0 + 28.0 /num_spins, 1e-10);
}