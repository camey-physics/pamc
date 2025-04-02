#include "models/TestModel.hpp"
#include <stdexcept>

TestModel::TestModel(double state) : state_(state) {}

void TestModel::initializeState() {
    state_ = 0.0;
}
void TestModel::initializeState(double state) {
    state_ = state;
}
void TestModel::copyStateFrom(const Model& other) {
    const TestModel* otherTestModel = dynamic_cast<const TestModel*>(&other);
    if (!otherTestModel) {
        throw std::invalid_argument("Cannot copy state from a non-TestModel instance.");
    }
    this->state_ = otherTestModel->state_;
}

double TestModel::calcEnergy() const {
    return state_ * state_;  // Example: Energy as squared state
}

void TestModel::updateSweep(int numSweeps) {
    state_ += numSweeps * 0.1;  // Example: Arbitrary update rule
}