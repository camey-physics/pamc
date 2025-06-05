#include <gtest/gtest.h>
#include <gsl/gsl_rng.h>
#include <vector>

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
    gsl_rng_set(rng_, 42);

    pop_size = 10;
    population = std::make_unique<Population<IsingModel>>(
        pop_size, gsl_rng_mt19937, *shared_data);
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