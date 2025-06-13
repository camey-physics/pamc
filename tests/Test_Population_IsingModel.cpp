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
    L = 4;
    num_spins = L * L * L;
    num_neighbors = 6;

    neighbor_table = initializeNeighborTable3D(L);
    bond_table.resize(num_spins * num_neighbors, -1.0);  // ferromagnetic bonds

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

TEST_F(PopulationIsingModelTest, RunEquilibrationSweep) {
  std::vector<std::vector<int>> initialState(pop_size);
  
  for (int i = 0; i < pop_size; ++i) {
    initialState[i] = population->getState(i);
  }

  population->equilibrate(1, 0.0, IsingModel::UpdateMethod::metropolis, true);

  for (int i = 0; i < pop_size; ++i) {
    EXPECT_NE(initialState[i], population->getState(i));
  }
}