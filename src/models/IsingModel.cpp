#include "IsingModel.hpp"

#include <gsl/gsl_rng.h>

#include <cassert>
#include <cmath>
#include <stdexcept>

IsingModel::IsingModel(const SharedModelData<IsingModel>& shared_data)
    : num_spins_(shared_data.num_spins),
      system_size_(shared_data.system_size),
      neighbor_table_(shared_data.neighbor_table),
      bond_table_(shared_data.bond_table) {
  spins_ = new int8_t[num_spins_];

  initializeState();  // Randomize the spins or set them as per the class design
}