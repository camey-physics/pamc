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

#include "Model.hpp"
#include "SharedModelData.hpp"

template <typename ModelType>
class Population {
 public:
  Population(unsigned int pop_size, const gsl_rng_type* T,
             const SharedModelData<ModelType>& shared_data);
  ~Population();
  template <typename... Args>
  void equilibrate(double beta, int num_sweeps,
                   typename ModelType::UpdateMethod method,
                   Args&&... extra_args);
  void resample(double new_beta);
  // Keep beta schedule simple for now.
  double suggestNextBeta() { return beta_ + 0.05; }
  double measureEnergy(bool force = false);

  // Not enforced to be in derived models via Model.hpp, but required for
  // Population. void measureObservable(typename ModelType::Observable);

  // population_[i].getState() is duck typed and not enforced by Model.hpp.
  auto getState(int i) const { return population_[i].getState(); }
  int const getBeta() const { return beta_; }
  double const getDeltaBetaF() const { return delta_betaF; }
  int const getPopSize() const { return pop_size_; }

  // Returns a const reference to the population of models for direct
  // interaction when needed. Not intended to be used in normal circumstances;
  // for unit testing and debugging.
  // const std::vector<ModelType>& getModels() const { return population_; }
  std::vector<ModelType>& getModels() { return population_; }

 private:
  double beta_ = 0.0;
  double delta_betaF_ = 0.0;
  unsigned int pop_size_ = 0;
  unsigned int nom_pop_size_ = 0;
  unsigned int max_pop_size_ = 0;
  std::vector<ModelType> population_;
  std::vector<double> energies_;
  std::vector<double> weights_;
  bool energies_current_ = false;
  const SharedModelData<ModelType>& shared_data_;

  gsl_rng* r_ = nullptr;

  // Helper functions
  void resizePopulationStorage(unsigned int new_size);
  inline int stochastic_round(double tau, gsl_rng* r) {
    int base = static_cast<int>(std::floor(tau));
    double prob = tau - base;
    return (gsl_rng_uniform(r) < prob) ? base + 1 : base;
  }
};

template <typename ModelType>
Population<ModelType>::Population(unsigned int pop_size, const gsl_rng_type* T,
                                  const SharedModelData<ModelType>& shared_data)
    : beta_(0.0),
      pop_size_(pop_size),
      nom_pop_size_(pop_size),
      shared_data_(shared_data),
      r_(gsl_rng_alloc(T)) {
  max_pop_size_ =
      static_cast<int>(nom_pop_size_ + 10 * std::sqrt(nom_pop_size_));
  resizePopulationStorage(pop_size);
  for (unsigned int i = 0; i < pop_size_; ++i) {
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
  for (unsigned int i = 0; i < pop_size_; ++i) {
    population_[i].updateSweep(num_sweeps, method, beta,
                               std::forward<Args>(extra_args)...);
  }
  energies_current_ = false;
}

template <typename ModelType>
void Population<ModelType>::resample(double new_beta) {
  double delta_beta = new_beta - beta_;
  double avg_energy = measureEnergy();

  // Apply energy shift to stabilize exponentials:
  // We compute weights ∝ exp(-Δβ (E_i - ⟨E⟩)) to avoid underflow,
  // and account for the shift in Q and betaF update.
  for (unsigned int i = 0; i < pop_size_; ++i) {
    weights_[i] = std::exp(-delta_beta * (energies_[i] - avg_energy));
  }
  double QR = std::accumulate(weights_.begin(), weights_.end(), 0.0);
  // Requires the additional part to cancel out the shifted energy
  delta_betaF_ -= std::log(QR / pop_size_) + delta_beta * avg_energy;
  // Now normalize weights (equal to tau_i).
  // The shifted energy in Q and in weights cancel each other.
  for (unsigned int i = 0; i < pop_size_; ++i) {
    weights_[i] = nom_pop_size_ * weights_[i] / QR;
  }
  // Implement actual resampling here.

  unsigned int total_new = 0;
  for (unsigned int i = 0; i < pop_size_; ++i) {
    weights_[i] = stochastic_round(weights_[i], r_);
    total_new += static_cast<int>(weights_[i]);
  }

  // Make space if needed
  if (total_new > population_.capacity()) {
    resizePopulationStorage(total_new);
  }

  unsigned int copy_to = 0;
  // Find the first zero-copy slot
  while (copy_to < pop_size_ && weights_[copy_to] > 0) {
    ++copy_to;
  }

  unsigned int copy_from = 0;
  while (copy_from < pop_size_) {
    double& count = weights_[copy_from];

    if (count <= 1) {
      ++copy_from;
      continue;
    }

    // Perform a copy
    population_[copy_to].copyStateFrom(population_[copy_from]);
    energies_[copy_to] = energies_[copy_from];
    --count;

    // Move copy_to to next available zero slot
    do {
      ++copy_to;
    } while (copy_to < pop_size_ && weights_[copy_to] > 0);

    if (copy_to >= pop_size_)
      break;
  }

  if (total_new < pop_size_) {
    copy_to = 0;
    copy_from = pop_size_ - 1;

    while (copy_to < copy_from) {
      // Find next zero-copy slot
      while (copy_to < copy_from && weights_[copy_to] > 0) ++copy_to;

      // Find next copyable replica from the end
      while (copy_to < copy_from && weights_[copy_from] == 0) --copy_from;

      if (copy_to < copy_from) {
        population_[copy_to].copyStateFrom(population_[copy_from]);
        energies_[copy_to] = energies_[copy_from];
        ++copy_to;
        --copy_from;
      }
    }

    population_.resize(total_new);
    energies_.resize(total_new);
    weights_.resize(total_new);
    pop_size_ = total_new;
  }
}

template <typename ModelType>
double Population<ModelType>::measureEnergy(bool force) {
  if (!energies_current_ || force) {
    for (unsigned int i = 0; i < pop_size_; ++i) {
      energies_[i] = population_[i].measureEnergy();
    }
    energies_current_ = true;
  }
  double total_energy =
      std::accumulate(energies_.begin(), energies_.end(), 0.0);
  return total_energy / pop_size_;
}

template <typename ModelType>
void Population<ModelType>::resizePopulationStorage(unsigned int new_size) {
  if (new_size > max_pop_size_) {
    throw std::runtime_error("Exceeded maximum allowed population size.");
  } else if (new_size > population_.capacity()) {
    int reserve_size = static_cast<int>(
        new_size + 5 * std::sqrt(static_cast<double>(new_size)));
    population_.reserve(reserve_size);
    energies_.reserve(reserve_size);
    weights_.reserve(reserve_size);
  }
  population_.resize(new_size);
  energies_.resize(new_size);
  weights_.resize(new_size);
}

#endif  // POPULATION_HPP