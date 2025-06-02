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

#include "Model.hpp"
#include "SharedModelData.hpp"

#include <gsl/gsl_rng.h>
#include <vector>

template <typename ModelType>
class Population {
 public:
  Population(int pop_size, const gsl_rng_type* T, const SharedModelData<ModelType>& shared_data);
  ~Population();
  template <typename... Args>
  void equilibrate(double beta, int num_sweeps, typename ModelType::UpdateMethod method, Args&&... extra_args);
  void resample(double new_beta);
  double suggestNextBeta();
  double measureEnergy();

  // Not enforced to be in derived models via Model.hpp, but required for Population.
  // void measureObservable(typename ModelType::Observable);

  // population_[i].getState() is duck typed and not enforced by Model.hpp.
  auto getState(int i) const { return population_[i].getState(); }

  // Returns a const reference to the population of models for direct interaction when needed. Not intended to be used in normal circumstances; for unit testing and debugging.
  const std::vector<ModelType>& getModels() const { return population_; }

 private:
  double beta_ = 0.0;
  int pop_size_ = 0;
  std::vector<ModelType> population_;
  const SharedModelData<ModelType>& shared_data_;

  gsl_rng* r_ = nullptr;
};

template <typename ModelType>
Population<ModelType>::Population(int pop_size, const gsl_rng_type* T, const SharedModelData<ModelType>& shared_data)
    : beta_(0.0),
      pop_size_(pop_size),
      population_(pop_size),
      r_(gsl_rng_alloc(T)),
      shared_data_(shared_data) {
  for (int i = 0; i < pop_size_; ++i) {
    population_[i].initializeState(r_, shared_data);
  }
}

template <typename ModelType>
Population<ModelType>::~Population() {
  gsl_rng_free(r_);
}

// Use a variadic template in order to pass additional arguments to the default arguments defined in Model.hpp
template <typename ModelType>
template <typename... Args>
void Population<ModelType>::equilibrate(double beta, int num_sweeps, typename ModelType::UpdateMethod method, Args&&... extra_args) {
  beta_ = beta;
  for (int i = 0; i < pop_size_; ++i) {
    population_[i].updateSweep(num_sweeps, method, beta, std::forward<Args>(extra_args)...);
  }
}

#endif // POPULATION_HPP