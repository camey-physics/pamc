#include <gsl/gsl_rng.h>
#include <gtest/gtest.h>
#include <math.h>

#include <vector>

#include "SharedModelData.hpp"
#include "models/Ising3DHelpers.hpp"
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

TEST(IsingModelTest, CopyState) {
  const int L = 10;
  const int num_spins = L * L * L;
  const int num_neighbors = 6;
  double J = 1;

  std::vector<int> neighbor_table = initializeNeighborTable3D(L);
  std::vector<double> bond_table(num_spins * num_neighbors, J);

  SharedModelData<IsingModel> data(L, num_spins, num_neighbors,
                                   neighbor_table.data(), bond_table.data());
  IsingModel model(data);
  IsingModel model2(data);
  gsl_rng* r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r, 42);
  model.initializeState(r);
  model2.copyStateFrom(model);

  for (int i = 0; i < num_spins; ++i) {
    EXPECT_EQ(model.getSpin(i), model2.getSpin(i));
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
  std::vector<int> expectedNeighbors = {(L - 1) * L * L, L * L, (L - 1) * L, L,
                                        L - 1,           1};
  for (int n = 0; n < num_neighbors; ++n) {
    EXPECT_EQ(table[spin * num_neighbors + n], expectedNeighbors[n]);
  }

  expectedNeighbors = {0, 32, 28, 20, 19, 17};
  spin = L * L;
  for (int n = 0; n < num_neighbors; ++n) {
    EXPECT_EQ(table[spin * num_neighbors + n], expectedNeighbors[n]);
  }
}

TEST(IsingModelTest, CalcEnergy) {
  const int L = 5;
  const int num_spins = L * L * L;
  const int num_neighbors = 6;
  double J = 1.0;

  std::vector<int> neighbor_table = initializeNeighborTable3D(L);
  std::vector<double> bond_table(num_spins * num_neighbors, J);

  SharedModelData<IsingModel> data(L, num_spins, num_neighbors,
                                   neighbor_table.data(), bond_table.data());

  // Calculate energy of all spins up configuration
  IsingModel model(data);
  EXPECT_NEAR(model.calcEnergy(), -3.0, 1e-10);
  // Calculate energy with one spin down
  int ind = index3D(1, 0, 0, L);
  model.setSpin(ind, -1);
  EXPECT_NEAR(model.calcEnergy(), -3.0 + 12.0 / num_spins, 1e-10);
  // Calculate energy with neighboring spin down
  ind = index3D(0, 0, 0, L);
  model.setSpin(ind, -1);
  EXPECT_NEAR(model.calcEnergy(), -3.0 + 20.0 / num_spins, 1e-10);
  // Calculate energy with second neighboring spin down
  ind = index3D(0, L - 1, 0, L);
  model.setSpin(ind, -1);
  EXPECT_NEAR(model.calcEnergy(), -3.0 + 28.0 / num_spins, 1e-10);
}

// Check that the metropolis algorithm obtains expected high temperature results
TEST(IsingModelTest, MetropolisSweep) {
  const int L = 10;
  const int num_spins = L * L * L;
  const int num_neighbors = 6;
  double beta = 0.1;
  int num_samples = 100;
  double J = 1;

  std::vector<int> neighbor_table = initializeNeighborTable3D(L);
  std::vector<double> bond_table(num_spins * num_neighbors, J);

  SharedModelData<IsingModel> data(L, num_spins, num_neighbors,
                                   neighbor_table.data(), bond_table.data());
  IsingModel model(data);
  gsl_rng* r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r, 42);
  model.initializeState(r);
  model.updateSweep(1000, beta, r, IsingModel::UpdateMethod::metropolis, true);

  double avg_energy = 0.0, avg_mag = 0.0;
  for (int i = 0; i < num_samples; ++i) {
    model.updateSweep(1000, beta, r, IsingModel::UpdateMethod::metropolis, true);
    avg_energy += model.calcEnergy();
    avg_mag += model.calcMagnetization();
  }
  avg_energy /= num_samples;
  avg_mag /= num_samples;

  EXPECT_NEAR(avg_energy, -3 * J * tanh(beta * J), 5e-2);
  EXPECT_NEAR(avg_mag, 0.0, 5e-2);
}

// Check that the heat bath algorithm obtains expected high temperature results
TEST(IsingModelTest, HeatBathSweep) {
  const int L = 10;
  const int num_spins = L * L * L;
  const int num_neighbors = 6;
  double beta = 0.1;
  int num_samples = 100;
  double J = 1;

  std::vector<int> neighbor_table = initializeNeighborTable3D(L);
  std::vector<double> bond_table(num_spins * num_neighbors, J);

  SharedModelData<IsingModel> data(L, num_spins, num_neighbors,
                                   neighbor_table.data(), bond_table.data());
  IsingModel model(data);
  gsl_rng* r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r, 42);
  model.initializeState(r);
  model.updateSweep(1000, beta, r, IsingModel::UpdateMethod::heat_bath, true);

  double avg_energy = 0.0, avg_mag = 0.0;
  for (int i = 0; i < num_samples; ++i) {
    model.updateSweep(1000, beta, r, IsingModel::UpdateMethod::metropolis, true);
    avg_energy += model.calcEnergy();
    avg_mag += model.calcMagnetization();
  }
  avg_energy /= num_samples;
  avg_mag /= num_samples;

  EXPECT_NEAR(avg_energy, -3 * J * tanh(beta * J), 5e-2);
  EXPECT_NEAR(avg_mag, 0.0, 5e-2);
}

// Check that the Wolff algorithm obtains expected low temperature results
TEST(IsingModelTest, WolffSweep) {
  const int L = 10;
  const int num_spins = L * L * L;
  const int num_neighbors = 6;
  double beta = 10;
  int num_samples = 100;
  double J = 1;

  std::vector<int> neighbor_table = initializeNeighborTable3D(L);
  std::vector<double> bond_table(num_spins * num_neighbors, J);

  SharedModelData<IsingModel> data(L, num_spins, num_neighbors,
                                   neighbor_table.data(), bond_table.data());
  IsingModel model(data);
  gsl_rng* r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r, 42);
  model.initializeState(r);
  model.updateSweep(20, beta, r, IsingModel::UpdateMethod::wolff, false);

  double avg_energy = 0.0, avg_mag = 0.0;
  for (int i = 0; i < num_samples; ++i) {
    model.updateSweep(10, beta, r, IsingModel::UpdateMethod::wolff, false);
    avg_energy += model.calcEnergy();
    avg_mag += abs(model.calcMagnetization());
  }
  avg_energy /= num_samples;
  avg_mag /= num_samples;

  EXPECT_NEAR(avg_energy, -3, 5e-2);
  EXPECT_NEAR(avg_mag, 1.0, 5e-2);
}