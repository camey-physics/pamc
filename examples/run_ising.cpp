#include <iostream>
#include <cstdlib>
#include <vector>
#include <gsl/gsl_rng.h>

#include "Population.hpp"
#include "models/IsingModel.hpp"
#include "SharedModelData.hpp"
#include "models/Ising3DHelpers.hpp"
#include "Genealogy.hpp"

int main(int argc, char* argv[]) {
    if (argc != 6) {
        std::cerr << "Usage: " << argv[0] << " <L> <pop_size> <culling_frac> <beta_max> <seed>" << std::endl;
        return 1;
    }

    // Parse arguments
    int L = std::atoi(argv[1]);
    int pop_size = std::atoi(argv[2]);
    double culling_frac = std::atof(argv[3]);
    double beta_min = 0.0;
    double beta_max = std::atof(argv[4]);
    unsigned long int seed = static_cast<unsigned long int>(std::stoul(argv[5]));

    // Prepare shared data
    int num_spins = L * L * L;
    int num_neighbors = 6;
    std::vector<int> neighbor_table = initializeNeighborTable3D(L);
    std::vector<double> bond_table(num_spins * num_neighbors, 1.0);

    SharedModelData<IsingModel> shared_data(L, num_spins, num_neighbors, neighbor_table.data(), bond_table.data());

    // Create population
    Population<IsingModel> population(pop_size, gsl_rng_mt19937, shared_data, seed);

    // Annealing loop
    double beta = beta_min;
    int step = 0;
    while (beta <= beta_max) {
        population.equilibrate(10, beta, IsingModel::UpdateMethod::metropolis, true);
        double E = population.measureEnergy(); 

        double M_sum = 0.0;
        double M2_sum = 0.0;
        double M4_sum = 0.0;
        for (const auto& model : population.getModels()) {
            double m = model.measureMagnetization() / num_spins;  
            M_sum += m;
            M2_sum += m * m;
            M4_sum += m * m * m * m;
        }
        int current_pop_size = population.getPopSize();
        double M_avg = M_sum / current_pop_size;
        double M2_avg = M2_sum / current_pop_size;
        double M4_avg = M4_sum / current_pop_size;

        double binder = 1.0 - M4_avg / (3.0 * M2_avg * M2_avg);

        GenealogyStatistics stats = population.computeGenealogyStatistics();

        std::cout << step << " "
        << beta << " " 
        << E / num_spins << " "
        << M_avg << " "
        << binder << " "
        << stats.rho_t << " "
        << stats.rho_s << std::endl;
        if (beta == beta_max) break;
        beta = population.suggestNextBeta(beta, culling_frac);
        if (beta > beta_max) beta = beta_max;
        population.resample(beta);
        step++;
    }

    return 0;
}
