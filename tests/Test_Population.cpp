#include <gtest/gtest.h>
#include <gsl/gsl_rng.h>

#include "Population.hpp"
#include "models/TestModel.hpp"
#include "SharedModelData.hpp"

class PopulationTestModelTest : public ::testing::Test {
 protected:
  void SetUp() override {
    gsl_rng_env_setup();
    rng_ = gsl_rng_alloc(gsl_rng_taus);

    pop_size = 5;
    shared_data = SharedModelData<TestModel>{};
    pop = std::make_unique<Population<TestModel>>(pop_size, gsl_rng_taus, shared_data);
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