// Population.cpp
#include "Population.hpp"
#include "TestModel.hpp"
#include <cmath>

template <typename ModelType>
Population<ModelType>::Population(int numModels, int modelSize) {
    int tmp = sqrt(numModels) + numModels;
    population_.reserve(numModels); // Preallocate memory for models
    for (int i = 0; i < numModels; ++i) {
        ModelType model;
        model.initialize(modelSize); // Initialize each model with size
        population_.push_back(std::move(model)); // Move model into the vector
    }
}

template <typename ModelType>
Population<ModelType>::~Population() {
    // No need to manually delete anything, as the vector will handle this.
}

template <typename ModelType>
void Population<ModelType>::update() {
    // Sweep over the entire population and update each model
    for (auto& model : population_) {
        model.updateSweep();  // Update each individual model
    }
}

template <typename ModelType>
double Population<ModelType>::measureEnergy() {
    // Calculate the average energy of the population
    double totalEnergy = 0.0;
    for (auto& model : population_) {
        totalEnergy += model.computeEnergy();
    }
    return totalEnergy / population_.size();
}

// template <typename ModelType>
// double Population<ModelType>::measureMagnetization() {
//     // Calculate the average magnetization of the population
//     double totalMagnetization = 0.0;
//     for (auto& model : population_) {
//         totalMagnetization += model.computeMagnetization();
//     }
//     return totalMagnetization / population_.size();
// }

// Explicit template instantiations for specific model types
template class Population<TestModel>;
