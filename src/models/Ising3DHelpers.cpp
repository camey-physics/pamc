#include "models/Ising3DHelpers.hpp"

std::vector<int> initializeNeighborTable3D(int system_size) {
  int L = system_size;
  int num_spins = L * L * L;
  std::vector<int> neighbor_table(num_spins * 6);

  for (int i = 0; i < L; ++i) {
    for (int j = 0; j < L; ++j) {
      for (int k = 0; k < L; ++k) {
        int ind = index3D(i, j, k, L);
        neighbor_table[ind * 6 + 0] = index3D(i - 1, j, k, L);
        neighbor_table[ind * 6 + 1] = index3D(i + 1, j, k, L);
        neighbor_table[ind * 6 + 2] = index3D(i, j - 1, k, L);
        neighbor_table[ind * 6 + 3] = index3D(i, j + 1, k, L);
        neighbor_table[ind * 6 + 4] = index3D(i, j, k - 1, L);
        neighbor_table[ind * 6 + 5] = index3D(i, j, k + 1, L);
      }
    }
  }

  return neighbor_table;
}
// Future arbitrary dimension support
// // Compute flattened index from D-dimensional coordinates
// int indexND(const std::vector<int>& coords, int L) {
//     int index = 0;
//     int stride = 1;
//     for (int d = coords.size() - 1; d >= 0; --d) {
//         int coord = (coords[d] % L + L) % L;  // periodic
//         index += coord * stride;
//         stride *= L;
//     }
//     return index;
// }

// // Compute D-dimensional neighbor table
// std::vector<int> makeNeighborTable(int L, int D) {
//     const int num_sites = std::pow(L, D);
//     const int num_neighbors = 2 * D;
//     std::vector<int> neighbor_table(num_sites * num_neighbors);

//     for (int flat = 0; flat < num_sites; ++flat) {
//         // Convert flat index to coordinates
//         std::vector<int> coords(D);
//         int temp = flat;
//         for (int d = D - 1; d >= 0; --d) {
//             coords[d] = temp % L;
//             temp /= L;
//         }

//         // Add ±1 neighbors in each dimension
//         for (int d = 0; d < D; ++d) {
//             std::vector<int> plus = coords;
//             std::vector<int> minus = coords;
//             plus[d] = (coords[d] + 1) % L;
//             minus[d] = (coords[d] - 1 + L) % L;

//             neighbor_table[flat * num_neighbors + 2 * d]     = indexND(minus,
//             L); neighbor_table[flat * num_neighbors + 2 * d + 1] =
//             indexND(plus, L);
//         }
//     }

//     return neighbor_table;
// }
