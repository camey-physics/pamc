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
  void initializeState() override;
  void copyStateFrom(const Model& other) override;

  // enum class UpdateMethod {};
  double calcEnergy() const override;
  // Required by the Model.hpp interface, not needed for IsingModel.
  void updateSweep(int num_sweeps) override {
    throw std::logic_error(
        "IsingModel::updateSweep(num_sweeps) is not supported — use "
        "updateSweep(num_sweeps, beta) instead.");
  }
  void updateSweep(int num_sweeps, double beta);

 private:
  int num_spins_;
  int system_size_;
  const int* neighbor_table_;
  const double* bond_table_;
  int8_t* spins_;
};

#endif  // ISING_MODEL_HPP
