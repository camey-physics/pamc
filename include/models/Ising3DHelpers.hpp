#ifndef ISING_3D_HELPERS_HPP
#define ISING_3D_HELPERS_HPP

#include <vector>

// Periodic boundary helper
inline int mod(int i, int system_size) {
  return (i % system_size + system_size) % system_size;
}

// Flatten 3D (i,j,k) index to 1D for a cubic system
inline int index3D(int i, int j, int k, int system_size) {
  return mod(i, system_size) * system_size * system_size +
         mod(j, system_size) * system_size + mod(k, system_size);
}

// Initialize 3D nearest-neighbor table for a cubic lattice with periodic
// boundaries. Returns a flat vector of size (system_size^3 × 6), where each
// spin has 6 neighbors.
std::vector<int> initializeNeighborTable3D(int system_size);

#endif  // ISING_3D_HELPERS_HPP