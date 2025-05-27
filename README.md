# Population Annealing Monte Carlo (PAMC)

This project is a clean, modular C++ implementation of the Population Annealing Monte Carlo (PAMC) algorithm, with an initial focus on spin glass systems in the canonical ensemble. It is designed to be extensible, testable, and performant, with shared disorder structures, parallelism, and long-term flexibility for scientific and industrial applications.

## Project Goals

- Provide a flexible simulation framework for population annealing
- Support clean separation of per-replica state and shared model data
- Apply modern C++ techniques for safety, clarity, and performance
- Enable future extensions to distributed (MPI) or GPU backends

## Core Features and Roadmap

- [x] Clean `Model` base class
- [x] Initial `IsingModel` implementation with shared `ModelData`
- [x] `Population` class to manage replicas, annealing, and resampling
- [ ] Resampling mechanism (multinomial to start)
- [ ] Observable tracking with statistical error estimators
- [ ] Temperature schedule generation (linear, geometric, constant culling fraction)
- [ ] OpenMP-based parallel update sweeps

## Code Structure (In Progress)

```
include/ # Public headers
src/ # Core implementation
models/ # Model definitions (e.g. IsingModel)
tests/ # Unit tests (GoogleTest)
```

## License

This project is licensed under the GNU General Public License v3.0 (GPLv3). Future versions may adopt a more permissive license to support broader adoption in academic and industrial settings.