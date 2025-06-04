#include <gtest/gtest.h>
#include <gsl/gsl_rng.h>
#include <iostream>

#include "Population.hpp"
#include "models/TestModel.hpp"
#include "SharedModelData.hpp"

class PopulationTestModelTest : public ::testing::Test {
 protected:
  void SetUp() override {
    gsl_rng_env_setup();
    rng_ = gsl_rng_alloc(gsl_rng_mt19937);

    pop_size = 10;
    shared_data = SharedModelData<TestModel>{};
    pop = std::make_unique<Population<TestModel>>(pop_size, gsl_rng_mt19937, shared_data);
  }

  void TearDown() override {
    gsl_rng_free(rng_);
  }

  int pop_size;
  gsl_rng* rng_ = nullptr;
  SharedModelData<TestModel> shared_data;
  std::unique_ptr<Population<TestModel>> pop;
};

TEST_F(PopulationTestModelTest, EquilibrateWithFakeMidYieldsEnergyInExpectedRange) {
  pop->equilibrate(1.0, 3, TestModel::UpdateMethod::FAKE_MID, rng_);

  for (int i = 0; i < pop_size; ++i) {
    double e = pop->getState(i);
    EXPECT_GE(e, 1.0);
    EXPECT_LT(e, 2.0);
  }

  double avg_e = pop->measureEnergy();
  EXPECT_GE(avg_e, 1.0);
  EXPECT_LT(avg_e, 2.0);

  const auto& models = pop->getModels();
  for (const auto& m : models) {
    EXPECT_TRUE(m.state_initialized);
    EXPECT_EQ(m.updates_called_, 3);
  }
}

TEST_F(PopulationTestModelTest, ResampleWithTwoEnergiesZeroTemperatureLimit) {
  std::vector<TestModel>& models = pop->getModels();
  for (int i = 0; i < 5; ++i) {
    models[i].setState(0.1);
  }
  for (int i = 5; i < 10; ++i) {
    models[i].setState(1);
  }
  pop->resample(10);
  EXPECT_EQ(pop->getPopSize(), 10);
  for (int i = 0; i < 10; ++i) {
    EXPECT_EQ(models[i].getState(), 0.1);
  }
}

// Here we explicitly test that tau_i = avg(n_i), where tau is the normalized resampling weight 
// and n is the stochastically determined number of copies. If this test fails, look at how large
// the discrepancy is and see if increasing num_trials fixes the problem.
TEST_F(PopulationTestModelTest, ResampleWithMultipleEnergiesCalculateTau) {
  // std::vector<double> expected_tau = {1.505449880326551, 1.362187382697220, 1.232558114240915, 1.115264701669020, 
  //                                     1.009133233084841, 0.913101509078768, 0.826208411879570, 0.747584286164701, 
  //                                     0.676442235257524, 0.612070245600891};
  std::vector<double> expected_tau = {1.20927569, 1.0941979, 0.9900712, 0.89585347, 0.81060174};

  // Explicitly set the resampling RNG simple reproducibility
  gsl_rng* r_local = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r_local, 7777);
  pop_size = 5;
  std::vector<int> n(pop_size, 0);
  int num_trials = 20000;
  int avg_pop_size = 0;

  // Reset the population for each trial, then resample and measure number of copies (n)
  for (int trial = 0; trial < num_trials; ++trial) {
    pop = std::make_unique<Population<TestModel>>(pop_size, gsl_rng_taus, shared_data);
    pop->setRngSeed(1001+trial);
    std::vector<TestModel>& models = pop->getModels();
    for (int i = 0; i < pop_size; ++i) {
      models[i].setState(1 + 0.1 *i);
      models[i].setFamily(i);
    }
    pop->resample(1, r_local);
    models = pop->getModels();
    for (int i = 0; i < pop->getPopSize(); ++i) {
      EXPECT_EQ(pop->getPopSize(), models.size());
      int ind = models[i].getFamily();
      n[ind] += 1;
    }
  }

  // Test that tau_i = avg(n_i) (averaged over resamples)
  for (int i = 0; i < pop_size; ++i) {
    EXPECT_NEAR(static_cast<double>(n[i]) /num_trials, expected_tau[i], 1e-2);
  }
  gsl_rng_free(r_local);
}

// TEST_F(PopulationTestModelTest, ResampleSeveralTimesMeasureFamily) {
//   // Explicitly set the resampling RNG simple reproducibility
//   gsl_rng* r_local = gsl_rng_alloc(gsl_rng_mt19937);
//   gsl_rng_set(r_local, 7777);
//   pop_size = 10;
//   int num_trials = 20000;
//   int avg_pop_size = 0;

//   // Reset the population for each trial, then resample and measure number of copies (n)
//   for (int trial = 0; trial < num_trials; ++trial) {
//     pop = std::make_unique<Population<TestModel>>(pop_size, gsl_rng_taus, shared_data);
//     pop->setRngSeed(1001+trial);
//     std::vector<TestModel>& models = pop->getModels();
//     for (int i = 0; i < pop_size; ++i) {
//       models[i].setState(1 + 0.1 *i);
//       models[i].setFamily(i);
//     }
//     pop->resample(1, r_local);
//     models = pop->getModels();
//     for (int i = 0; i < pop->getPopSize(); ++i) {
//       EXPECT_EQ(pop->getPopSize(), models.size());
//       int ind = models[i].getFamily();
//       n[ind] += 1;
//     }
//   }

//   // Test that tau_i = avg(n_i) (averaged over resamples)
//   for (int i = 0; i < pop_size; ++i) {
//     EXPECT_NEAR(static_cast<double>(n[i]) /num_trials, expected_tau[i], 1e-2);
//   }
//   gsl_rng_free(r_local);
// }