#ifndef GENEOLOGY_HPP
#define GENEOLOGY_HPP

struct GenealogyStatistics {
    explicit GenealogyStatistics(int init_size)
        : initial_pop_size(init_size) {}

    const int initial_pop_size;
    double rho_t = 0.0;
    double rho_s = 0.0;
    double culling_frac_target = 0.0;
    double culling_frac_actual = 0.0;
    int num_unique_families = 0;
    int max_family_size = 0;
};

#endif