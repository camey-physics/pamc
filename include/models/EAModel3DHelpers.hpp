#ifndef EA_3D_HELPERS_HPP
#define EA_3D_HELPERS_HPP

#include <vector>
#include <string>

// Load a 3D lattice neighbor table from file (same format as Ising3DHelpers)
std::vector<int> loadNeighborTable(const std::string& filename, int num_spins, int num_neighbors);

// Load bond values from file (each line = site, columns = bonds)
std::vector<double> loadBondTable(const std::string& filename, int num_spins, int num_neighbors);

#endif  // EA_MODEL_3D_HELPERS_HPP