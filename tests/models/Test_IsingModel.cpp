#include <gsl/gsl_rng.h>
#include <gtest/gtest.h>
#include <cmath>

#include <vector>

#include "SharedModelData.hpp"
#include "models/Ising3DHelpers.hpp"
#include "models/IsingModel.hpp"

class TestIsingModel : public ::testing::Test {
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

TEST_F(TestIsingModel, Construct) {
  IsingModel model(shared_data);

  for (int i = 0; i < num_spins; ++i) {
    EXPECT_EQ(model.getSpin(i), 1);
  }
}

TEST_F(TestIsingModel, CopyState) {
  IsingModel model(shared_data);
  IsingModel model2(shared_data);
  gsl_rng* r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r, 42);
  model.initializeState(r);
  model2.copyStateFrom(model);

  for (int i = 0; i < num_spins; ++i) {
    EXPECT_EQ(model.getSpin(i), model2.getSpin(i));
  }
  gsl_rng_free(r);
}

TEST_F(TestIsingModel, GetState) {
  IsingModel model(shared_data);

  gsl_rng* r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r, 42);
  model.initializeState(r);
  std::vector<int> state = model.getState();

  EXPECT_EQ(num_spins, state.size());
  
  for (int i = 0; i < num_spins; ++i) {
    EXPECT_EQ(model.getSpin(i), state[i]);
  }
  gsl_rng_free(r);
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

TEST_F(TestIsingModel, MeasureEnergy) {
  // Measure energy of all spins up configuration
  IsingModel model(shared_data);
  EXPECT_NEAR(model.measureEnergy(), -3.0 *num_spins, 1e-10);
  // Measure energy with one spin down
  int ind = index3D(1, 0, 0, L);
  model.setSpin(ind, -1);
  EXPECT_NEAR(model.measureEnergy(), -3.0 *num_spins + 12.0, 1e-10);
  // Measure energy with neighboring spin down
  ind = index3D(0, 0, 0, L);
  model.setSpin(ind, -1);
  EXPECT_NEAR(model.measureEnergy(), -3.0 *num_spins + 20.0, 1e-10);
  // Measure energy with second neighboring spin down
  ind = index3D(0, L - 1, 0, L);
  model.setSpin(ind, -1);
  EXPECT_NEAR(model.measureEnergy(), -3.0 *num_spins + 28.0, 1e-10);
}

// Check that the metropolis algorithm obtains expected high temperature results
TEST_F(TestIsingModel, MetropolisSweep) {
  double beta = 0.1;
  int num_samples = 100;

  IsingModel model(shared_data);
  gsl_rng* r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r, 42);
  model.initializeState(r);
  model.updateSweep(1000, beta, r, IsingModel::UpdateMethod::metropolis, true);

  double avg_energy = 0.0, avg_mag = 0.0;
  for (int i = 0; i < num_samples; ++i) {
    model.updateSweep(1000, beta, r, IsingModel::UpdateMethod::metropolis,
                      true);
    avg_energy += model.measureEnergy();
    avg_mag += model.measureMagnetization();
  }
  avg_energy /= num_samples *num_spins;
  avg_mag /= num_samples *num_spins;

  EXPECT_NEAR(avg_energy, -3 * J * tanh(beta * J), 5e-2);
  EXPECT_NEAR(avg_mag, 0.0, 5e-2);

  gsl_rng_free(r);
}

// Check that the heat bath algorithm obtains expected high temperature results
TEST_F(TestIsingModel, HeatBathSweep) {
  double beta = 0.1;
  int num_samples = 100;

  IsingModel model(shared_data);
  gsl_rng* r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r, 42);
  model.initializeState(r);
  model.updateSweep(1000, beta, r, IsingModel::UpdateMethod::heat_bath, true);

  double avg_energy = 0.0, avg_mag = 0.0;
  for (int i = 0; i < num_samples; ++i) {
    model.updateSweep(1000, beta, r, IsingModel::UpdateMethod::metropolis,
                      true);
    avg_energy += model.measureEnergy();
    avg_mag += model.measureMagnetization();
  }
  avg_energy /= num_samples *num_spins;
  avg_mag /= num_samples *num_spins;

  EXPECT_NEAR(avg_energy, -3 * J * tanh(beta * J), 5e-2);
  EXPECT_NEAR(avg_mag, 0.0, 5e-2);
  gsl_rng_free(r);
}

// Check that the Wolff algorithm obtains expected low temperature results
TEST_F(TestIsingModel, WolffSweep) {
  double beta = 10;
  int num_samples = 100;

  IsingModel model(shared_data);
  gsl_rng* r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r, 42);
  model.initializeState(r);
  model.updateSweep(20, beta, r, IsingModel::UpdateMethod::wolff, false);

  double avg_energy = 0.0, avg_mag = 0.0;
  for (int i = 0; i < num_samples; ++i) {
    model.updateSweep(10, beta, r, IsingModel::UpdateMethod::wolff, false);
    avg_energy += model.measureEnergy();
    avg_mag += abs(model.measureMagnetization());
  }
  avg_energy /= num_samples *num_spins;
  avg_mag /= num_samples *num_spins;

  EXPECT_NEAR(avg_energy, -3, 5e-2);
  EXPECT_NEAR(avg_mag, 1.0, 5e-2);
  gsl_rng_free(r);
}