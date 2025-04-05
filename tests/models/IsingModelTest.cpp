#include "models/IsingModel.hpp"
#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <tuple>

class IsingModelTest : public ::testing::Test {
protected:
    int L = 4;
    double beta = 1.0;
    int seed = 5000;
    double J = 1.0;

    IsingModel createModel(int l = 4, double b = 1.0, double j = 1.0, int s = 5000) {
        return IsingModel(l, b, j, s);
    }
};

TEST_F(IsingModelTest, InitializesSpinLattice) {
    IsingModel model = createModel();
    for (int i = 0; i < L; ++i)
        for (int j = 0; j < L; ++j)
            for (int k = 0; k < L; ++k)
                EXPECT_TRUE(model.getSpin(i, j, k) == +1 || model.getSpin(i, j, k) == -1);
}

TEST_F(IsingModelTest, PeriodicBoundaryConditions) {
    IsingModel model = createModel();
    model.setSpin(L-1, L-1, L-1, -1);
    EXPECT_EQ(model.getSpin(-1, -1, -1), model.getSpin(L-1, L-1, L-1));
}

TEST_F(IsingModelTest, NeighborTable) {
    IsingModel model = createModel();
    std::vector<int> expectedNeighbors = {(L-1)*L*L, L*L, (L-1)*L, L, L-1, 1};
    EXPECT_EQ(model.getNeighbors(0), expectedNeighbors);

    expectedNeighbors = {0, 32, 28, 20, 19, 17};
    EXPECT_EQ(model.getNeighbors(L * L), expectedNeighbors);
}

TEST_F(IsingModelTest, TestParams) {
    IsingModel model = createModel();
    auto [L_out, J_out, beta_out, seed_out] = model.getParams();
    EXPECT_EQ(L_out, L);
    EXPECT_EQ(J_out, J);
    EXPECT_EQ(beta_out, beta);
    EXPECT_EQ(seed_out, seed);
}

TEST_F(IsingModelTest, CalcEnergy) {
    int L5 = 5;
    IsingModel model(L5, 1, 1, 5000);
    model.initializeAllUp();
    EXPECT_NEAR(model.calcEnergy(), -3.0, 1e-10);

    model.setSpin(1, 0, 0, -1);
    EXPECT_NEAR(model.calcEnergy(), -3.0 + 12.0 / (L5 * L5 * L5), 1e-10);

    model.setSpin(0, 0, 0, -1);
    EXPECT_NEAR(model.calcEnergy(), -3.0 + 20.0 / (L5 * L5 * L5), 1e-10);

    model.setSpin(0, L5 - 1, 0, -1);
    EXPECT_NEAR(model.calcEnergy(), -3.0 + 28.0 / (L5 * L5 * L5), 1e-10);
}

TEST_F(IsingModelTest, MetropolisSweep) {
    int L10 = 10;
    int numSamples = 100;
    double smallBeta = 0.1;

    IsingModel model(L10, smallBeta, J, seed);
    double avg_energy = 0.0, avg_mag = 0.0;

    model.updateSweep(1000, IsingModel::UpdateMethod::metropolis, true);

    for (int i = 0; i < numSamples; ++i) {
        model.updateSweep(100, IsingModel::UpdateMethod::metropolis, true);
        avg_energy += model.calcEnergy();
        avg_mag += model.calcMagnetization();
    }

    avg_energy /= numSamples;
    avg_mag /= numSamples;

    EXPECT_NEAR(avg_energy, -3 *J *tanh(smallBeta *J), 5e-2);
    EXPECT_NEAR(avg_mag, 0.0, 5e-2);
}

TEST_F(IsingModelTest, HeatBathSweep) {
    int L10 = 10;
    int numSamples = 100;
    double smallBeta = 0.1;

    IsingModel model(L10, smallBeta, J, seed);
    model.initializeAllUp();
    double avg_energy = 0.0, avg_mag = 0.0;

    model.updateSweep(1000, IsingModel::UpdateMethod::heatBath, true);

    for (int i = 0; i < numSamples; ++i) {
        model.updateSweep(100, IsingModel::UpdateMethod::heatBath, true);
        avg_energy += model.calcEnergy();
        avg_mag += model.calcMagnetization();
    }

    avg_energy /= numSamples;
    avg_mag /= numSamples;

    EXPECT_NEAR(avg_energy, -3 *J *tanh(smallBeta *J), 5e-2);
    EXPECT_NEAR(avg_mag, 0.0, 5e-2);
}

TEST_F(IsingModelTest, LowTempWolff) {
    int L10 = 10;
    int numSamples = 100;
    double smallBeta = 10;

    IsingModel model(L10, smallBeta, J, seed);
    model.initializeAllUp();
    double avg_energy = 0.0, avg_mag = 0.0;

    model.updateSweep(20, IsingModel::UpdateMethod::wolff);

    for (int i = 0; i < numSamples; ++i) {
        model.updateSweep(10, IsingModel::UpdateMethod::wolff);
        avg_energy += model.calcEnergy();
        avg_mag += model.calcMagnetization();
    }

    avg_energy /= numSamples;
    avg_mag /= numSamples;

    EXPECT_NEAR(avg_energy, -3, 5e-2);
    EXPECT_NEAR(abs(avg_mag), 1.0, 5e-2);
}

TEST_F(IsingModelTest, CopyState) {
    int L1 = 10;
    int L2 = 10;
    int seed1 = 497235;

    IsingModel model1(L1, 0.1, J, seed1);
    IsingModel model2(L2, beta, J, seed);

    model1.updateSweep(100, IsingModel::UpdateMethod::heatBath, true);
    model2.copyStateFrom(model1);

    EXPECT_EQ(model1.calcEnergy(), model2.calcEnergy());

    auto [Lcheck1, Jcheck1, betacheck1, seedcheck1] = model1.getParams();
    auto [Lcheck2, Jcheck2, betacheck2, seedcheck2] = model2.getParams();

    EXPECT_EQ(Lcheck1, Lcheck2);
    EXPECT_EQ(Jcheck1, Jcheck2);
    EXPECT_EQ(betacheck1, betacheck2);
    EXPECT_EQ(seedcheck1, seed1);
    EXPECT_EQ(seedcheck2, seed);
}


int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
