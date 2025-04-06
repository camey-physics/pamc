#ifndef ISING_MODEL_H
#define ISING_MODEL_H

#include "models/Model.hpp"
#include <vector>
#include <gsl/gsl_rng.h>
#include <memory>
#include <tuple>

class IsingModel : public Model {
public:
    explicit IsingModel(int L, double beta, double J, int seed);
    IsingModel(int L, double beta, double J, int seed, int* externalSpins);
    ~IsingModel();
    void initializeState() override; // Randomizes spins
    void initializeAllUp();
    void copyStateFrom(const Model& other) override;
    enum class UpdateMethod { metropolis, heatBath, wolff };

    // Observer methods
    int getSpin(int i, int j, int k) const;
    std::vector<int> getNeighbors(int index) const;
    double getBeta() const;
    auto getParams() const {
        return std::make_tuple(L_, J_, beta_, seed_);
    }
    int* getSpinPointer() const;

    // Modifier methods
    void setSpin(int i, int j, int k, int val);
    void setBeta(double beta);

    // Computation methods
    double calcEnergy() const override;
    double calcMagnetization() const;

    void updateSweep(int numSweeps) override {
        updateSweep(numSweeps, UpdateMethod::metropolis, false);
    }

    void updateSweep(int numSweeps, UpdateMethod method, bool sequential=0);

private:
    int L_;
    double J_;
    double beta_;
    int seed_;
    int* spins_;
    bool ownsSpins_ = true;
    std::vector<int> NT_;

    // Monte Carlo update methods
    void metropolis(int i);
    void heatBath(int i);
    int wolff();

    inline int index(int i, int j, int k) const;
    inline int mod(int i) const;
    void initializeNT();
    gsl_rng *r;
    int calcLocalH(int i) const;

};

#endif // ISING_MODEL_H
