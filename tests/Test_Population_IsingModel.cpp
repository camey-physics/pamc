#include <gtest/gtest.h>
#include <gsl/gsl_rng.h>
#include <vector>
#include <cmath>

#include "Population.hpp"
#include "models/IsingModel.hpp"
#include "SharedModelData.hpp"
#include "models/Ising3DHelpers.hpp"

class PopulationIsingModelTest : public ::testing::Test {
 protected:
  void SetUp() override {
    L = 5;
    num_spins = L * L * L;
    num_neighbors = 6;
    J = 1.0;

    neighbor_table = initializeNeighborTable3D(L);
    bond_table.resize(num_spins * num_neighbors, J);

    shared_data = new SharedModelData<IsingModel>(L, num_spins, num_neighbors, neighbor_table.data(), bond_table.data());

    rng_ = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng_, 56);

    pop_size = 1000;
    population = std::make_unique<Population<IsingModel>>(
        pop_size, gsl_rng_mt19937, *shared_data, 6416);
  }

  void TearDown() override {
    gsl_rng_free(rng_);
    delete shared_data;
  }

  int L;
  int num_spins;
  int num_neighbors;
  double J;
  int pop_size;
  gsl_rng* rng_;
  std::vector<int> neighbor_table;
  std::vector<double> bond_table;
  SharedModelData<IsingModel>* shared_data;
  std::unique_ptr<Population<IsingModel>> population;
};

class LargePopulationIsingModelTest : public ::testing::Test {
 protected:
  void SetUp() override {
    L = 5;
    num_spins = L * L * L;
    num_neighbors = 6;
    J = 1.0;

    neighbor_table = initializeNeighborTable3D(L);
    bond_table.resize(num_spins * num_neighbors, J);  // ferromagnetic bonds

    shared_data = new SharedModelData<IsingModel>(L, num_spins, num_neighbors, neighbor_table.data(), bond_table.data());

    rng_ = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng_, 56);

    pop_size = 10;
    population = std::make_unique<Population<IsingModel>>(
        pop_size, gsl_rng_mt19937, *shared_data, 6416);
  }

  void TearDown() override {
    gsl_rng_free(rng_);
    delete shared_data;
  }

  int L;
  int num_spins;
  int num_neighbors;
  double J;
  int pop_size;
  gsl_rng* rng_;
  std::vector<int> neighbor_table;
  std::vector<double> bond_table;
  SharedModelData<IsingModel>* shared_data;
  std::unique_ptr<Population<IsingModel>> population;
};

TEST_F(PopulationIsingModelTest, PopulationInitializesCorrectly) {
  std::vector<IsingModel>& models = population->getModels();
  EXPECT_EQ(models.size(), pop_size);
}

TEST_F(PopulationIsingModelTest, GetState) {
  std::vector<IsingModel>& models = population->getModels();
  for (int i = 0; i < pop_size; ++i) {
    auto state = population->getState(i);
    EXPECT_EQ(state, models[i].getState());
  }
}

// Test that the initial configuration (spins randomly oriented) has an energy within two
// standard deviations that of the expected value of E = 0.
TEST_F(PopulationIsingModelTest, MeasureInitialEnergy) {
  EXPECT_NEAR(population->measureEnergy(), 0.0, sqrt(3.0 /num_spins));
}

TEST_F(PopulationIsingModelTest, RunEquilibrationSweepAtBetaZero) {
  std::vector<std::vector<int>> initialState(pop_size);
  
  for (int i = 0; i < pop_size; ++i) {
    initialState[i] = population->getState(i);
  }

  population->equilibrate(10, 0.0, IsingModel::UpdateMethod::metropolis, false);

  for (int i = 0; i < pop_size; ++i) {
    EXPECT_NE(initialState[i], population->getState(i));
  }
  // The energy per spin should still be close to zero
  EXPECT_NEAR(population->measureEnergy(), 0.0, sqrt(3.0 /num_spins));
}

// Test annealing with metropolis and heat bath for large beta.
TEST_F(LargePopulationIsingModelTest, AnnealWithMetropolisToLowBeta) {
  double beta = 0.05;
  while (beta <= 0.15) {
    population->equilibrate(200, beta, IsingModel::UpdateMethod::metropolis, true);
    EXPECT_NEAR(population->measureEnergy(), -3 * J * tanh(beta * J), 5e-2);
    beta += 0.05;
    population->resample(beta);
  }
}

TEST_F(LargePopulationIsingModelTest, AnnealWithHeatBathToLowBeta) {
  double beta = 0.05;
  while (beta <= 0.15) {
    population->equilibrate(200, beta, IsingModel::UpdateMethod::heat_bath, false);
    EXPECT_NEAR(population->measureEnergy(), -3 * J * tanh(beta * J), 5e-2);
    beta += 0.05;
    population->resample(beta);
  }
}

TEST_F(LargePopulationIsingModelTest, EquilibrateWithWolffAtHighBeta) {
  double beta = 10;
  population->equilibrate(100, beta, IsingModel::UpdateMethod::wolff, false);
  EXPECT_NEAR(population->measureEnergy(), -3, 5e-2);
}