# Population Annealing Monte Carlo (PAMC)

Flexible Population Annealing Monte Carlo (PAMC) codebase for general spin models and optimization problems. Currently supports fully functioning 3D Ising simulations, with planned extensions to spin glasses, QUBO optimization, and high-performance computing.

---

## Project Status

This repository is under active development. Core functionality is implemented for population annealing on the 3D Ising model, with adaptive temperature schedules and verified basic observables. The code is structured for modular extension and intended for both research and pedagogical purposes.

---

## Core Features

### Algorithm

- [x] Core model infrastructure (`Model` interface)
- [x] `IsingModel` class with `Metropolis`, `heat_bath`, and `Wolff` updates
- [x] `Population` class for managing replicas, annealing, and resampling
- [x] Resampling mechanism (multinomial resampling)
- [x] Adaptive temperature schedule (`Population::suggestNextBeta()` using energy variance)
- [x] Verified Binder cumulant crossing (3D Ising) for Population/Observable infrastructure
- [x] Genealogical observables (`rho_t`, `rho_s`, max family size, etc)

### Infrastructure & Postprocessing

- [x] Unit tests (GoogleTest) and continuous integration
- [ ] Generic observable interface (`measureObservable()`)
- [ ] Data output and I/O framework (`DataWriter` class or equivalent)
- [ ] Python postprocessing and analysis scripts for examples and tests
- [ ] OpenMP-based parallel update sweeps

---

## Directory Layout

- `include/` — Public headers
  - `Model.hpp` — abstract model interface
  - `Population.hpp` — population annealing engine
  - `SharedModelData.hpp` — shared model parameters (neighbor tables, bond tables)
  - `models/` — model-specific headers (e.g. `IsingModel.hpp`, `TestModel.hpp`)
- `src/` — Model implementations (e.g. `models/IsingModel.cpp`)
- `examples/` — Standalone simulation drivers (e.g. `run_ising.cpp`)
- `tests/` — Unit tests (GoogleTest)
- `validation/` — Python scripts for validating simulation output and generating analysis plots (e.g., Binder cumulant crossing)

---

## Planned Models

- Fully connected and lattice-based spin glasses
- Flexible shared data structures for generalized coupling matrices

---

## Examples

Example simulations can be found in `examples/`. For instance, the `run_ising.cpp` program performs adaptive annealing runs and outputs data suitable for Binder cumulant crossing analysis.

---

## Build Instructions

CMake is used to configure and build the project.

### Debug build (with address sanitizer and extra warnings):

```bash
cmake -DCMAKE_BUILD_TYPE=Debug -B build-debug
cmake --build build-debug
```

### Release build:

```bash
cmake -DCMAKE_BUILD_TYPE=Release -B build-release
cmake --build build-release
```
