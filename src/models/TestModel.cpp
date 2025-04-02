#include "models/TestModel.hpp"

TestModel::TestModel(double state) : state_(state) {}

// virtual ~Model() = default;
// virtual void initializeState() = 0;
// virtual void copyStateFrom(const Model& other) = 0;
// enum class UpdateMethod { };
// virtual double calcEnergy() const = 0;
// virtual void updateSweep(int numSweeps) = 0;