#ifndef TEST_MODEL_HPP
#define TEST_MODEL_HPP

#include <gsl/gsl_rng.h>

#include <stdexcept>
#include "Model.hpp"
#include "SharedModelData.hpp"

class TestModel {
 public:
  TestModel(const SharedModelData<TestModel>&) {}
  enum class UpdateMethod { FAKE_LOW, FAKE_MID, FAKE_HIGH };

  void initializeState(gsl_rng* r) {
    state_initialized = true;
    energy_ = gsl_rng_uniform(r);  // [0,1)
  }

  // void updateSweep(int num_sweeps, UpdateMethod method, double /*beta*/, gsl_rng* r) {
  void updateSweep(int num_sweeps, double /*beta*/, gsl_rng* r,
                             UpdateMethod method, bool /*sequential*/) {
    for (int i = 0; i < num_sweeps; ++i) {
      updates_called_ += 1;
    }

    switch (method) {
      case UpdateMethod::FAKE_LOW:
        energy_ = gsl_rng_uniform(r);             // [0,1)
        break;
      case UpdateMethod::FAKE_MID:
        energy_ = 1.0 + gsl_rng_uniform(r);       // [1,2)
        break;
      case UpdateMethod::FAKE_HIGH:
        energy_ = 2.0 + gsl_rng_uniform(r);       // [2,3)
        break;
      default:
        throw std::invalid_argument("Unknown update method in TestModel");
    }
  }

  double measureEnergy() const {
    return energy_;
  }

  double getState() const {
    return energy_;
  }

  void copyStateFrom(const TestModel& other) {
    energy_ = other.energy_;
    updates_called_ = other.updates_called_;
    state_initialized = other.state_initialized;
    family_ = other.family_;
    parent_ = other.parent_;
  }

  void setState(double energy) { energy_ = energy; }
  void setFamily(int family) {
    if (family_ != -1) {
      throw std::logic_error("family_ already set");
    }
    family_ = family;
  }
  void setParent(int parent) { parent_ = parent; }

  int getFamily() const { return family_; }
  int getParent() const { return parent_; }

  bool state_initialized = false;
  int updates_called_ = 0;

 private:
  double energy_ = 0.0;
  int family_ = -1;
  int parent_ = -1;
};

#endif  // TEST_MODEL_HPP