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
  IsingModel(const SharedModelData<IsingModel>& shared_data, int* external_spins);
  ~IsingModel();

  // Static helper functions related to external memory
  static constexpr bool supports_pool = true;
  using storage_type = int;
  static std::size_t elementsPerReplica(const SharedModelData<IsingModel>& shared_data) { return shared_data.num_spins; }

  void initializeState(gsl_rng* r) override;
  void copyStateFrom(const Model& other) override;

  // IsingModel specific enumerated classes
  enum class UpdateMethod { metropolis, heat_bath, wolff }; 
  enum class Observable { energy, magnetization };

  // IsingModel observable methods
  double measureEnergy() const override;
  double measureMagnetization() const;

  // Monte Carlo sweep methods.
  // By default uses metropolis updates on randomly selected spins
  void updateSweep(int num_sweeps, double beta, gsl_rng* r) override {
    updateSweep(num_sweeps, beta, r, UpdateMethod::metropolis, false);
  }
  void updateSweep(int num_sweeps, double beta, gsl_rng* r, UpdateMethod method,
                   bool sequential = false);

  const std::vector<int> getState() const { return std::vector<int>(spins_, spins_ + num_spins_); }

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

  // Helper methods for unit testing IsingModel class
  void setSpin(int i, int val);
  int getSpin(int i) const;

  int* getSpins_() const { return spins_; }

 private:
  // Shared model data, immutable
  const int num_spins_;
  const int num_neighbors_;
  const int system_size_;
  const int* neighbor_table_;
  const double* bond_table_;
  bool owns_spins_ = true;
  int family_ = -1;
  int parent_ = -1;
  int* spins_;

  // Monte Carlo update methods
  void metropolis(gsl_rng* r, double beta, int i);
  void heatBath(gsl_rng* r, double beta, int i);
  int wolff(gsl_rng* r, double beta);

  // Helper functions
};

#endif  // ISING_MODEL_HPP
