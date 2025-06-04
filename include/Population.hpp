#ifndef POPULATION_HPP
#define POPULATION_HPP

// Population<ModelType> assumes that ModelType defines the following methods:
//
//   void initializeState(gsl_rng*, const SharedModelData<ModelType>&);
//   void updateSweep(int, UpdateMethod, double, gsl_rng*, bool [, ...]);
//   auto getState() const;         // return type may vary
//
// These methods are duck-typed: they are not enforced via Model.hpp,
// but are required for Population to compile and function correctly.

#include <gsl/gsl_rng.h>

#include <cmath>
#include <numeric>
#include <vector>
#include <iostream>

#include "Model.hpp"
#include "SharedModelData.hpp"

template <typename ModelType>
class Population {
 public:
  Population(int pop_size, const gsl_rng_type* T,
             const SharedModelData<ModelType>& shared_data);
  ~Population();
  template <typename... Args>
  void equilibrate(double beta, int num_sweeps,
                   typename ModelType::UpdateMethod method,
                   Args&&... extra_args);
  void resample(double new_beta, gsl_rng* r_override = nullptr);
  // Keep beta schedule simple for now.
  double suggestNextBeta() { return beta_ + 0.05; }
  double measureEnergy(bool force = false);

  // Not enforced to be in derived models via Model.hpp, but required for
  // Population. void measureObservable(typename ModelType::Observable);

  // population_[i].getState() is duck typed and not enforced by Model.hpp.
  auto getState(int i) const { return population_[i].getState(); }
  int getBeta() const { return beta_; }
  double getDeltaBetaF() const { return delta_betaF_; }
  int getPopSize() const { return pop_size_; }

  void setRngSeed(unsigned long int s) { gsl_rng_set(r_, s); }
  void setNomPopSize(int i) { nom_pop_size_ = i; }

  // Returns a const reference to the population of models for direct
  // interaction when needed. Not intended to be used in normal circumstances;
  // for unit testing and debugging.
  std::vector<ModelType>& getModels() { return population_; }

 private:
  double beta_ = 0.0;
  double delta_betaF_ = 0.0;
  int pop_size_ = 0;
  int nom_pop_size_ = 0;
  int max_pop_size_ = 0;
  std::vector<ModelType> population_;
  std::vector<double> energies_;
  std::vector<double> weights_;
  std::vector<int> copy_counts_;
  bool energies_current_ = false;
  const SharedModelData<ModelType>& shared_data_;

  gsl_rng* r_ = nullptr;

  // Helper functions
  void resizePopulationStorage(int new_size);
  inline int stochastic_round(double tau, gsl_rng* r) {
    int floor = static_cast<int>(std::floor(tau));
    double prob = tau - floor;
    double rng = gsl_rng_uniform(r);
    return (rng < prob) ? floor + 1 : floor;
  }
  void computeWeights(double new_beta, double avg_energy, double& QR);
  void computeCopyCounts(int& total_new, gsl_rng* r_local);
  void forwardCopy(int old_pop_size, int new_pop_size);
  void backfillHoles(int old_pop_size);
};

template <typename ModelType>
Population<ModelType>::Population(int pop_size, const gsl_rng_type* T,
                                  const SharedModelData<ModelType>& shared_data)
    : beta_(0.0),
      pop_size_(pop_size),
      nom_pop_size_(pop_size),
      shared_data_(shared_data),
      r_(gsl_rng_alloc(T)) {
  max_pop_size_ =
      static_cast<int>(nom_pop_size_ + 10 * std::sqrt(nom_pop_size_));
  resizePopulationStorage(pop_size);
  for (int i = 0; i < pop_size_; ++i) {
    population_[i].initializeState(r_, shared_data);
  }
}

template <typename ModelType>
Population<ModelType>::~Population() {
  gsl_rng_free(r_);
}

// Use a variadic template in order to pass additional arguments to the default
// arguments defined in Model.hpp
template <typename ModelType>
template <typename... Args>
void Population<ModelType>::equilibrate(double beta, int num_sweeps,
                                        typename ModelType::UpdateMethod method,
                                        Args&&... extra_args) {
  beta_ = beta;
  for (int i = 0; i < pop_size_; ++i) {
    population_[i].updateSweep(num_sweeps, method, beta,
                               std::forward<Args>(extra_args)...);
  }
  energies_current_ = false;
}

// Resample to new_beta. Internally uses local copies of old_pop_size and new_pop_size
// to make logic clearer, since pop_size_ is updated indirectly by helpers.
template <typename ModelType>
void Population<ModelType>::resample(double new_beta, gsl_rng* r_override) {
  gsl_rng* r_local = r_override ? r_override : r_;
  double delta_beta = new_beta - beta_;
  double avg_energy = measureEnergy();
  double QR = 0.0;
  int old_pop_size = pop_size_;

  // This function calculates normalized weights (held in weights_) 
  // using a shifted energy for numerical stability. It also updates QR,
  // which is the normalization factor used in calculating delta (beta *F),
  // where the delta_beta *avg_energy term compensates for the energy shift
  computeWeights(new_beta, avg_energy, QR);
  delta_betaF_ -= std::log(QR / pop_size_) + delta_beta *avg_energy;

  // computeCopyCounts() updates both copy_counts_ and new_pop_size
  int new_pop_size = 0;
  computeCopyCounts(new_pop_size, r_local);

  if (new_pop_size >= old_pop_size) {
    resizePopulationStorage(new_pop_size);
    forwardCopy(old_pop_size, new_pop_size);
  }
  else {
    forwardCopy(old_pop_size, new_pop_size);
    backfillHoles(old_pop_size);
    resizePopulationStorage(new_pop_size);
  }
  pop_size_ = new_pop_size;

  int check_sum = std::accumulate(copy_counts_.begin(), copy_counts_.end(), 0);
  assert(check_sum == new_pop_size);
}

template <typename ModelType>
double Population<ModelType>::measureEnergy(bool force) {
  if (!energies_current_ || force) {
    for (int i = 0; i < pop_size_; ++i) {
      energies_[i] = population_[i].measureEnergy();
    }
    energies_current_ = true;
  }
  double total_energy =
      std::accumulate(energies_.begin(), energies_.end(), 0.0);
  return total_energy / pop_size_;
}

template <typename ModelType>
void Population<ModelType>::resizePopulationStorage(int new_size) {
  if (new_size > max_pop_size_) {
    throw std::runtime_error("Exceeded maximum allowed population size.");
  } else if (new_size > static_cast<int>(population_.capacity())) {
    int reserve_size = static_cast<int>(
        new_size + 5 * std::sqrt(static_cast<double>(new_size)));
    population_.reserve(reserve_size);
    energies_.reserve(reserve_size);
    weights_.reserve(reserve_size);
    copy_counts_.reserve(reserve_size);
  }
  population_.resize(new_size);
  energies_.resize(new_size);
  weights_.resize(new_size);
  copy_counts_.resize(new_size);
  for (int i = pop_size_; i < new_size; ++i) {
    energies_[i] = 0.0;
    weights_[i] = 0.0;
    copy_counts_[i] = 0;
  }
  pop_size_ = new_size;
}

// Computes QR (shifted) and the normalized weights tau (held in weights_)
template <typename ModelType>
void Population<ModelType>::computeWeights(double new_beta, double avg_energy, double& QR) {
  double delta_beta = new_beta - beta_;

  // Apply energy shift to stabilize exponentials:
  // We compute weights ∝ exp(-Δβ (E_i - ⟨E⟩)) to avoid underflow,
  // and account for the shift in the delta betaF update.
  for (int i = 0; i < pop_size_; ++i) {
    weights_[i] = std::exp(-delta_beta * (energies_[i] - avg_energy));
  }
  QR = std::accumulate(weights_.begin(), weights_.end(), 0.0);

  // Now normalize weights (equal to tau_i).
  // The shifted energy in Q and in weights cancel each other.
  for (int i = 0; i < pop_size_; ++i) {
    weights_[i] = nom_pop_size_ * weights_[i] / QR;
  }
}

template <typename ModelType>
void Population<ModelType>::computeCopyCounts(int& new_pop_size, gsl_rng* r_local) {
  new_pop_size = 0;
  for (int i = 0; i < pop_size_; ++i) {
    int n = stochastic_round(weights_[i], r_local);
    copy_counts_[i] = n;
    new_pop_size += n;
  }
}

template <typename ModelType>
void Population<ModelType>::forwardCopy(int old_pop_size, int new_pop_size) {
  int copy_from = 0;
  int copy_to = 0;

  // Search for the first index to copy to and to copy from
  while (copy_to < new_pop_size && copy_counts_[copy_to] > 0) ++copy_to;
  while (copy_from < old_pop_size && copy_counts_[copy_from] <= 1) ++copy_from;

  while (copy_from < old_pop_size && copy_to < new_pop_size) {
    population_[copy_to].copyStateFrom(population_[copy_from]);
    energies_[copy_to] = energies_[copy_from];

    --copy_counts_[copy_from];
    ++copy_counts_[copy_to];  // Optional; for debug or consistency checks

    while (copy_to < new_pop_size && copy_counts_[copy_to] > 0) ++copy_to;
    while (copy_from < old_pop_size && copy_counts_[copy_from] <= 1) ++copy_from;
  }
}

template <typename ModelType>
void Population<ModelType>::backfillHoles(int old_pop_size) {
  int copy_to = 0;
  int copy_from = old_pop_size - 1;

  // Move states from the back of the population to fill zero-copy slots in the front
  while (copy_to < copy_from) {
    // Advance to next hole in the front
    while (copy_to < copy_from && copy_counts_[copy_to] > 0) ++copy_to;

    // Find a model with at least 1 copy at the end
    while (copy_to < copy_from && copy_counts_[copy_from] == 0) --copy_from;

    if (copy_to < copy_from) {
      population_[copy_to].copyStateFrom(population_[copy_from]);
      energies_[copy_to] = energies_[copy_from];
      copy_counts_[copy_to] = 1;
      --copy_counts_[copy_from];
      ++copy_to;
      --copy_from;
    }
  }
}

#endif  // POPULATION_HPP