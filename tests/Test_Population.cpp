#include <gsl/gsl_rng.h>
#include <gtest/gtest.h>
#include <math.h>

#include "Population.hpp"
#include "models/TestModel.hpp"
#include "SharedModelData.hpp"

TEST(PopulationTest, EquilibrateCallsUpdateSweep) {
  SharedModelData<TestModel> dummy_data;  // could be empty
  Population<TestModel> pop(5, gsl_rng_taus, dummy_data);

  pop.equilibrate(1.0, 10, TestModel::UpdateMethod::FAKE);

  for (int i = 0; i < 5; ++i) {
    EXPECT_EQ(10,  pop.getState(i));
  }

  for (const auto& model : pop.getModels()) {
    EXPECT_TRUE(model.state_initialized);
    EXPECT_EQ(model.last_num_sweeps, 10);
    EXPECT_EQ(model.last_beta, 1.0);
  }
}
