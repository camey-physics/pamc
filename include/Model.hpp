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

  // Optional tag for model-specific observables (e.g. Energy, Magnetization,
  // Overlap).
  enum class Observables {};
};

#endif  // MODEL_HPP