#ifndef ISING_MODEL_HPP
#define ISING_MODEL_HPP

#include <gsl/gsl_rng.h>

#include <memory>
#include <stdexcept>
#include <vector>

#include "Model.hpp"
#include "SharedModelData.hpp"

class IsingModel : public Model {
 public:
  // IsingModel state related methods
  explicit IsingModel(const SharedModelData<IsingModel>& shared_data);
  ~IsingModel();
  void initializeState(gsl_rng* r) override;
  void copyStateFrom(const Model& other) override;

  // IsingModel specific enumerated classes
  enum class UpdateMethod { metropolis, heatBath, wolff };
  enum class Observables { energy, magnetization };

  // IsingModel observable methods
  double calcEnergy() const override;
  double calcMagnetization() const;


  void updateSweep(int num_sweeps, double beta, gsl_rng* r) override;

 private:
  // Shared model data, immutable
  const int num_spins_;
  const int num_neighbors_;
  const int system_size_;
  const int* neighbor_table_;
  const double* bond_table_;

  // Owned data
  int8_t* spins_;
};

#endif  // ISING_MODEL_HPP
