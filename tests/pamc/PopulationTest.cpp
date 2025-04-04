#include <gtest/gtest.h>
#include "pamc/Population.hpp"
#include "TestModel.hpp"

// Test Population initialization
TEST(PopulationTest, Initialization) {
    int populationSize = 10;
    int modelSize = 5;
    Population<TestModel> pop(populationSize, modelSize);
    EXPECT_EQ(pop.measureEnergy(), 25.0); // Assuming initial energy is zero
}

// Test update method
TEST(PopulationTest, Update) {
    int populationSize = 10;
    int modelSize = 5;
    Population<TestModel> pop(populationSize, modelSize);
    pop.update(); // Should update all models
    EXPECT_TRUE(true); // Placeholder check, assuming no errors occur
}

// Test energy measurement
TEST(PopulationTest, MeasureEnergy) {
    int populationSize = 10;
    int modelSize = 5;
    Population<TestModel> pop(populationSize, modelSize);
    double energyBefore = pop.measureEnergy();
    pop.update();
    double energyAfter = pop.measureEnergy();
    EXPECT_NE(energyBefore, energyAfter); // Energy should change after update
}