// Population.hpp
#ifndef POPULATION_HPP
#define POPULATION_HPP

#include <vector>
#include "Model.hpp"  // Base class Model

// Template definition for the Population class
template <typename ModelType>
class Population {
public:
    Population(int size, int modelSize = 0);
    ~Population();

    // Method to update the entire population (e.g., sweep over all models)
    void update();

    // Method to measure the average energy across the population
    double measureEnergy();

    // // Method to measure the average magnetization across the population
    // double measureMagnetization();

private:
    int size_;  // The size of the population
    std::vector<ModelType> population_;  // Vector to hold instances of the derived model type
};

#endif // POPULATION_HPP