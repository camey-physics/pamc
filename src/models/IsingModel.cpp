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
  for (int i = 0; i < num_spins_; ++i) {
    spins_[i] = 1;
  }
}

IsingModel::~IsingModel() { delete[] spins_; }

void IsingModel::initializeState(gsl_rng* r) {
  for (int i = 0; i < num_spins_; ++i) {
    int s = gsl_rng_uniform_int(r, 2) * 2 - 1;
    spins_[i] = s;
  }
}

void IsingModel::copyStateFrom(const Model& other) {
  const IsingModel& isingOther = static_cast<const IsingModel&>(other);
  assert(this->system_size_ == isingOther.system_size_ &&
         "System sizes must match!");
  assert(this->num_spins_ == isingOther.num_spins_ &&
         "Number of spins must match!");
  for (int i = 0; i < num_spins_; ++i) {
    this->spins_[i] = isingOther.spins_[i];
  }
}

double IsingModel::calcEnergy() const {
  double energy = 0.0;
  for (int i = 0; i < num_spins_; ++i) {
    for (int n = 0; n < 6; n += 2) {
      int j = neighbor_table_[i * 6 + n];
      energy += spins_[i] * spins_[j] * bond_table_[i * 6 + n];
    }
  }
  return energy / num_spins_;
}