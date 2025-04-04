#include <gtest/gtest.h>
#include "models/TestModel.hpp"
#include <cmath>

// Tolerance for floating-point comparisons.
const double kTol = 1e-5;

class TestModelTest : public ::testing::Test {
protected:
    TestModel* model;

    // Set up a TestModel with an initial state of 2.0.
    void SetUp() override {
        model = new TestModel(2.0);
    }

    void TearDown() override {
        delete model;
    }
};

// Test that the constructor initializes the state correctly.
// With an initial state of 2.0, calcEnergy() should return 4.0 (2.0^2).
TEST_F(TestModelTest, ConstructorAndCalcEnergy) {
    double energy = model->calcEnergy();
    EXPECT_NEAR(energy, 4.0, kTol);
}

// Test the default initializeState() method, which should set state to 0.0.
TEST_F(TestModelTest, InitializeStateDefault) {
    model->initializeState();
    // When state is 0.0, energy should be 0.0.
    double energy = model->calcEnergy();
    EXPECT_NEAR(energy, 0.0, kTol);
}

// Test the overloaded initializeState(double) method.
TEST_F(TestModelTest, InitializeStateWithParameter) {
    model->initializeState(3.0);
    // With state set to 3.0, energy should be 9.0.
    double energy = model->calcEnergy();
    EXPECT_NEAR(energy, 9.0, kTol);
}

// Test copyStateFrom by copying the state from another TestModel.
TEST_F(TestModelTest, CopyStateFrom) {
    // Create another TestModel with state 5.0.
    TestModel otherModel(5.0);
    model->copyStateFrom(otherModel);
    // After copying, the energy should be 25.0 (5.0^2).
    double energy = model->calcEnergy();
    EXPECT_NEAR(energy, 25.0, kTol);
}

// Test updateSweep to verify that the state updates correctly.
// Given that updateSweep(numSweeps) adds numSweeps*0.1 to the state,
// for an initial state of 2.0, updateSweep(10) should add 1.0, making the new state 3.0,
// so calcEnergy() should return 9.0.
TEST_F(TestModelTest, UpdateSweep) {
    model->updateSweep(10);
    double newEnergy = model->calcEnergy();
    EXPECT_NEAR(newEnergy, 9.0, kTol);
}
