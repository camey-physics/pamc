#include "models/IsingModel.hpp"

#include <gsl/gsl_rng.h>

#include <cassert>
#include <cmath>
#include <stdexcept>

IsingModel::IsingModel(const SharedModelData<IsingModel>& shared_data)
    : num_spins_(shared_data.num_spins),
      num_neighbors_(shared_data.num_neighbors),
      system_size_(shared_data.system_size),
      neighbor_table_(shared_data.neighbor_table),
      bond_table_(shared_data.bond_table) {
  assert(num_neighbors_ % 2 == 0 &&
         "Neighbor table must use even pairing (+/- directions)");
  spins_ = new int[num_spins_];
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
    this->family_ = isingOther.family_;
    this->parent_ = isingOther.parent_;
  }
}

double IsingModel::measureEnergy() const {
  double energy = 0.0;
  for (int i = 0; i < num_spins_; ++i) {
    // Skip every second neighbor to avoid double-counting bonds.
    // Assumes symmetric neighbor table with even num_neighbors_.
    for (int n = 0; n < num_neighbors_; n += 2) {
      int j = neighbor_table_[i * num_neighbors_ + n];
      energy -= spins_[i] * spins_[j] * bond_table_[i * num_neighbors_ + n];
    }
  }
  return energy / num_spins_;
}

double IsingModel::measureMagnetization() const {
  int mag = 0;
  for (int i = 0; i < num_spins_; ++i) {
    mag += spins_[i];
  }
  return static_cast<double>(mag) / num_spins_;
}

void IsingModel::updateSweep(int num_sweeps, double beta, gsl_rng* r,
                             UpdateMethod method, bool sequential) {
  void (IsingModel::*update_func)(gsl_rng*, double, int) = nullptr;
  switch (method) {
    case UpdateMethod::metropolis:
      update_func = &IsingModel::metropolis;
      break;
    case UpdateMethod::heat_bath:
      update_func = &IsingModel::heatBath;
      break;
    case UpdateMethod::wolff:
      if (sequential) {
        throw std::invalid_argument(
            "Wolff update cannot be used with sequential mode");
      }
      // wolff method handled separately below due to randomness in number of
      // flipped spins
      for (int sweep = 0; sweep < num_sweeps; ++sweep) {
        int num_flipped = 0;
        while (num_flipped < num_spins_) {
          num_flipped += wolff(r, beta);
        }
      }
      return;
    default:
      throw std::invalid_argument("Unknown update method!");
  }

  // For metropolis and heatBath, use the chosen update function pointer and
  // flip spins individually
  if (sequential) {
    for (int sweep = 0; sweep < num_sweeps; ++sweep) {
      for (int i = 0; i < num_spins_; ++i) {
        (this->*update_func)(r, beta, i);
      }
    }
  } else {
    for (int sweep = 0; sweep < num_sweeps; ++sweep) {
      for (int i = 0; i < num_spins_; ++i) {
        int s = gsl_rng_uniform_int(r, num_spins_);
        (this->*update_func)(r, beta, s);
      }
    }
  }
}

std::vector<int> IsingModel::getState() {
  std::vector<int> state(spins_, spins_ + num_spins_);
  return state;
}

void IsingModel::setSpin(int i, int val) {
  if (val != 1 && val != -1) {
    throw std::invalid_argument("Spin value must be +1 or -1");
  }
  spins_[i] = val;
}

int IsingModel::getSpin(int i) const {
  if (i < 0 || i >= num_spins_) {
    throw std::out_of_range("Index out of range");
  }
  return spins_[i];
}

void IsingModel::metropolis(gsl_rng* r, double beta, int i) {
  double delta_E = 0.0;
  for (int n = 0; n < num_neighbors_; ++n) {
    int j = neighbor_table_[i * num_neighbors_ + n];
    delta_E += spins_[j] * bond_table_[i * num_neighbors_ + n];
  }
  delta_E *= 2 * spins_[i];

  if (delta_E <= 0 || gsl_rng_uniform(r) < exp(-beta * delta_E)) {
    spins_[i] *= -1;
  }
}

void IsingModel::heatBath(gsl_rng* r, double beta, int i) {
  double local_h = 0.0;
  for (int n = 0; n < num_neighbors_; ++n) {
    int j = neighbor_table_[i * num_neighbors_ + n];
    local_h += spins_[j] * bond_table_[i * num_neighbors_ + n];
  }
  double probUp = 1 / (1 + exp(-2 * beta * local_h));
  if (gsl_rng_uniform(r) < probUp) {
    spins_[i] = 1;
  } else {
    spins_[i] = -1;
  }
}

int IsingModel::wolff(gsl_rng* r, double beta) {
  std::vector<bool> visited(num_spins_, false);
  std::vector<int> stack;
  // To keep track of the cluster size
  int clusterSize = 0;
  // This implementation of the Wolff algorithm is set up for Ising models that
  // have constant bond value.
  double J = bond_table_[0];

  // Pick a random starting spin
  int ind = gsl_rng_uniform_int(r, num_spins_);
  // Spin type of cluster
  int clusterSpin = spins_[ind];
  stack.push_back(ind);
  visited[ind] = true;

  // Probability of adding bond to cluster
  double P_add = 1 - exp(-2 * beta * J);

  while (!stack.empty()) {
    int i = stack.back();
    stack.pop_back();

    // Flip spin
    spins_[i] *= -1;
    clusterSize++;

    // Check neighbors
    for (int n = 0; n < num_neighbors_; ++n) {
      int j = neighbor_table_[i * num_neighbors_ + n];

      // If neighbor has the same spin and isn't visited, try adding to cluster
      if (!visited[j] && spins_[j] == clusterSpin &&
          gsl_rng_uniform(r) < P_add) {
        stack.push_back(j);
        visited[j] = true;
      }
    }
  }

  // Return the size of the cluster
  return clusterSize;
}