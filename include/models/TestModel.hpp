#ifndef TEST_MODEL_HPP
#define TEST_MODEL_HPP

#include <gsl/gsl_rng.h>

#include <memory>
#include <stdexcept>
#include <vector>

#include "Model.hpp"
#include "SharedModelData.hpp"

class TestModel {
 public:
  enum class UpdateMethod { FAKE };

  void initializeState(gsl_rng* r, const SharedModelData<TestModel>&) {
    state_initialized = true;
  }

  void updateSweep(int num_sweeps, UpdateMethod method, double beta) {
    last_num_sweeps = num_sweeps;
    last_beta = beta;
    for (int i = 0; i < num_sweeps; ++i) {
      updates_called += 1;
    }
  }

  double measureEnergy() const {
    return fixed_energy;
  }

  int getState() const { return updates_called; }

  // Optional: hooks to verify state
  bool state_initialized = false;
  int updates_called = 0;
  int last_num_sweeps = 0;
  double last_beta = 0.0;
  double fixed_energy = 42.0;  // default predictable result
};

#endif // TEST_MODEL_HPP