#ifndef MODEL_HPP
#define MODEL_HPP

#include <gsl/gsl_rng.h>

// #include "MemoryBlock.hpp"

// Abstract base class for all Monte Carlo models used in PAMC.
// Model is a polymorphic base for all models used in Population.
// SharedModelData<ModelType> is templated and not linked to Model directly;
// each derived model defines how it uses its own SharedModelData<ModelType>.

// Population<ModelType> assumes that ModelType defines the following methods:
//
//   void initializeState(gsl_rng*, const SharedModelData<ModelType>&);
//   void updateSweep(int, UpdateMethod, double, gsl_rng*, bool [, ...]);
//   auto getState() const;         // return type may vary
//
// These methods are duck-typed: they are not enforced via Model.hpp,
// but are required for Population to compile and function correctly.


class Model {
 public:
  virtual ~Model() = default;

  // Randomizes internal state using the provided RNG.
  // This base method must be overloaded in each derived class.
  // Models must also define a duck-typed overload that takes
  // SharedModelData<ModelType>&.
  virtual void initializeState(gsl_rng* r) = 0;

  // Copies the full model state from another instance.
  virtual void copyStateFrom(const Model& other) = 0;

  // Calculates and returns the current energy of the model.
  virtual double measureEnergy() const = 0;

  // Applies update sweeps with given beta and RNG.
  // Derived models may overload this to accept additional arguments.
  virtual void updateSweep(int num_sweeps, double beta, gsl_rng* r) = 0;

  // Optional tag for method-specific updates (e.g., Metropolis, Heat Bath).
  enum class UpdateMethod {};

  // Families can only be set once and is inherited via copyStateFrom
  void setFamily(int family) {
    if (family_ != -1) {
      throw std::logic_error("family_ already set");
    }
    family_ = family;
  }
  void setParent(int parent) { parent_ = parent; }

  int getFamily() const { return family_; }
  int getParent() const { return parent_; }

  // Duck-typed method: must be implemented in each model class and must take no
  // arguments.
  //   auto getState() const;
  // Used by Population to extract the state of a model instance.
  // The return type may vary across models.

  // Optional tag for model-specific observables (e.g., Energy, Magnetization).
  // enum class Observable {};

  // Generic observable measurement interface.
  // Still under design—will likely use duck typing:
  // Population expects a method measureObservable(...) where the argument
  // is ModelType::Observable (defined only in derived classes).
  private:
    int family_ = -1;
    int parent_ = -1;
};

#endif  // MODEL_HPP