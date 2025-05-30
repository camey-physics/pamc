#ifndef POPULATION_HPP
#define POPULATION_HPP

#include "Model.hpp"

#include <gsl/gsl_rng.h>
#include <vector>

template <typename ModelType>
class Population {
 public:
  Population(int pop_size, const gsl_rng_type* T);
  ~Population();
   
  void equilibrate(double beta, int numSweeps, typename ModelType::UpdateMethod method);
  void resample(double new_beta);
  double suggestNextBeta();
  double measureEnergy();

  // Not enforced to be in derived models via Model.hpp, but required for Population.
  // void measureObservable(typename ModelType::Observable);

 private:
  double beta_ = 0.0;
  int pop_size_ = 0;
  std::vector<ModelType> population_;

  gsl_rng* r_ = nullptr;
};

#endif // POPULATION_HPP