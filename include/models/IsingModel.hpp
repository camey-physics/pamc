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
  enum class UpdateMethod { metropolis, heat_bath, wolff };
  enum class Observable { energy, magnetization };

  // IsingModel observable methods
  double calcEnergy() const override;
  double calcMagnetization() const;

  // Monte Carlo sweep methods.
  // By default uses metropolis updates on randomly selected spins
  void updateSweep(int num_sweeps, double beta, gsl_rng* r) override {
    updateSweep(num_sweeps, beta, r, UpdateMethod::metropolis, false);
  }
  void updateSweep(int num_sweeps, double beta, gsl_rng* r, UpdateMethod method,
                   bool sequential = false);

  std::vector<int> getState();

  // Helper methods for unit testing IsingModel class
  void setSpin(int i, int val);
  int getSpin(int i) const;

 private:
  // Shared model data, immutable
  const int num_spins_;
  const int num_neighbors_;
  const int system_size_;
  const int* neighbor_table_;
  const double* bond_table_;

  // Owned data
  int* spins_;

  // Monte Carlo update methods
  void metropolis(gsl_rng* r, double beta, int i);
  void heatBath(gsl_rng* r, double beta, int i);
  int wolff(gsl_rng* r, double beta);

  // Helper functions
};

#endif  // ISING_MODEL_HPP
