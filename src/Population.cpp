#include <Population.hpp>

template <typename ModelType>
Population<ModelType>::Population(int pop_size, const gsl_rng_type* T)
    : beta_(0.0),
      pop_size_(pop_size),
      population_(pop_size),
      r_(gsl_rng_alloc(T)) {
  for (int i = 0; i < pop_size_; ++i) {
    population_[i].initializeState(r_);
  }
}

template <typename ModelType>
Population<ModelType>::~Population() {
  gsl_rng_free(r_);
}