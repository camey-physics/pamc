#ifndef SHARED_MODEL_DATA_HPP
#define SHARED_MODEL_DATA_HPP

// Primary template (unspecialized)
template <typename ModelT>
struct SharedModelData;

// Specialization for IsingModel
// The bond and neighbor tables are specified externally. The only constraint is
// that all spins must have the same number of neighbors.
// Consider changing neighbor_table and bond_table to std::span for bounds
// checking.
template<>
struct SharedModelData<class IsingModel> {
  const int system_size;
  const int num_spins;
  const int num_neighbors;
  const int* neighbor_table;
  const double* bond_table;
  SharedModelData(int system_size, int num_spins, int num_neighbors,
                  const int* neighbor_table, const double* bond_table)
      : system_size(system_size),
        num_spins(num_spins),
        num_neighbors(num_neighbors),
        neighbor_table(neighbor_table),
        bond_table(bond_table) {}
};

#endif