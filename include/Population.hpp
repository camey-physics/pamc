#ifndef POPULATION_HPP
#define POPULATION_HPP

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
      r_(gsl_rng_alloc(T)).
      shared_data_(shared_data) {
  for (int i = 0; i < pop_size_; ++i) {
    population_[i].initializeState(r_);
  }
}

template <typename ModelType>
Population<ModelType>::~Population() {
  gsl_rng_free(r_);
}

template <typename ModelType>
template <typename... Args>
void Population<ModelType>::equilibrate(double beta, int num_sweeps, typename ModelType::UpdateMethod method, Args&&... extra_args) {
  beta_ = beta;
  for (int i = 0; i < pop_size_; ++i) {
    population_[i].updateSweep(numSweeps, method, beta, std::forward<Args>(extra_args)...);
  }
}

#endif // POPULATION_HPP