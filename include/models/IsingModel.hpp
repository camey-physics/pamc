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
  explicit IsingModel(const SharedModelData<IsingModel>& shared_data);
  ~IsingModel();
  void initializeState(gsl_rng* r) override;
  void copyStateFrom(const Model& other) override;

  enum class UpdateMethod { metropolis, heatBath, wolff };
  double calcEnergy() const override;
  void updateSweep(int num_sweeps, double beta, gsl_rng* r) override;

 private:
  // Shared model data, immutable
  const int num_spins_;
  const int system_size_;
  const int* neighbor_table_;
  const double* bond_table_;
  int8_t* spins_;
};

#endif  // ISING_MODEL_HPP
