#ifndef SHARED_MODEL_DATA_HPP
#define SHARED_MODEL_DATA_HPP

// Primary template (unspecialized)
template <typename ModelT>
struct SharedModelData;

// Specialization for IsingModel
struct SharedModelData<class IsingModel> {
    int system_size;
    int num_spins;
    int num_neighbors;
    const int* neighbor_table;
    const double* bond_table;
};

#endif