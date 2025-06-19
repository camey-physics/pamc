#include <iostream>
#include <cstdlib>
#include <vector>
#include <string>
#include <gsl/gsl_rng.h>
#include <iomanip>

#include "Population.hpp"
#include "models/IsingModel.hpp"
#include "SharedModelData.hpp"
#include "models/EAModel3DHelpers.hpp"
#include "Genealogy.hpp"

int main(int argc, char* argv[]) {
    if (argc != 8) {
        std::cerr << "Usage: " << argv[0] 
                  << " <L> <pop_size> <culling_frac> <beta_max> <seed> <neighbor_table_path> <bond_table_path>" 
                  << std::endl;
        return 1;
    }

    // Parse arguments
    int L = std::atoi(argv[1]);
    int pop_size = std::atoi(argv[2]);
    double culling_frac = std::atof(argv[3]);
    double beta_max = std::atof(argv[4]);
    unsigned long int seed = static_cast<unsigned long int>(std::stoul(argv[5]));
    std::string neighbor_path = argv[6];
    std::string bond_path = argv[7];

    int num_spins = L * L * L;
    int num_neighbors = 6;

    std::vector<int> neighbor_table = loadNeighborTable(neighbor_path, num_spins, num_neighbors);
    std::vector<double> bond_table = loadBondTable(bond_path, num_spins, num_neighbors);

    SharedModelData<IsingModel> shared_data(L, num_spins, num_neighbors,
                                            neighbor_table.data(), bond_table.data());

    Population<IsingModel> population(pop_size, gsl_rng_mt19937, shared_data, seed);

    double beta = 0.0;
    int step = 0;
    while (beta <= beta_max) {
        population.equilibrate(10, beta, IsingModel::UpdateMethod::metropolis, true);
        double E = population.measureEnergy();
        double E_min = population.getMinEnergy();
        GenealogyStatistics stats = population.computeGenealogyStatistics();

        std::cout << std::fixed << std::setprecision(15)
          << step << " " << beta << " "
          << E << " "
          << E_min << " "
          << stats.rho_t << " "
          << stats.num_gs_families << std::endl;

        if (beta == beta_max) break;
        beta = population.suggestNextBeta(beta, culling_frac);
        if (beta > beta_max) beta = beta_max;
        population.resample(beta);
        step++;
    }

    return 0;
}
