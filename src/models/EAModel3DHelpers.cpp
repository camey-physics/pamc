#include <fstream>
#include <stdexcept>
#include <vector>
#include <string>

#include "models/EAModel3DHelpers.hpp"

std::vector<int> loadNeighborTable(const std::string& filename, int num_spins, int num_neighbors) {
    std::ifstream infile(filename);
    if (!infile) throw std::runtime_error("Failed to open neighbor table file.");

    std::vector<int> neighbors(num_spins * num_neighbors);
    for (int i = 0; i < num_spins; ++i) {
        for (int n = 0; n < num_neighbors; ++n) {
            infile >> neighbors[i * num_neighbors + n];
        }
    }
    return neighbors;
}

std::vector<double> loadBondTable(const std::string& filename, int num_spins, int num_neighbors) {
    std::ifstream infile(filename);
    if (!infile) throw std::runtime_error("Failed to open bond table file.");

    std::vector<double> bonds(num_spins * num_neighbors);
    for (int i = 0; i < num_spins; ++i) {
        for (int n = 0; n < num_neighbors; ++n) {
            infile >> bonds[i * num_neighbors + n];
        }
    }
    return bonds;
}