#ifndef MODEL_HPP
#define MODEL_HPP

#include <gsl/gsl_rng.h>

// Abstract base class for all Monte Carlo models used in PAMC.
class Model {
 public:
  virtual ~Model() = default;

  // Randomizes internal state using provided RNG.
  virtual void initializeState(gsl_rng* r) = 0;

  // Copies full model state from another instance.
  virtual void copyStateFrom(const Model& other) = 0;

  // Calculates and returns the current energy of the model.
  virtual double calcEnergy() const = 0;

  // Applies update sweeps with given beta and RNG.
  virtual void updateSweep(int num_sweeps, double beta, gsl_rng* r) = 0;

  // Optional tag for method-specific updates (e.g. Metropolis, Heat Bath).
  enum class UpdateMethod {};

  // Generic method to getState is duck typed because the return type
  // depends on the specific model. Should take no arguments, e.g.
  // std::vector<int> getState();

  // Optional tag for model-specific observables (e.g. Energy, Magnetization,
  // Overlap).
  // enum class Observable {};

  // Generic method to measure observables with generic input. Still need to
  // decide how to pass arbitrary observables through Population without it
  // needing to include specific model headers.
  // Currently considering duck typing so that Population expects a
  // measureObservable method with argument ModelType::observable, which is
  // only defined within the derived Model classes.
  // measureObservable(...)
};

#endif  // MODEL_HPP