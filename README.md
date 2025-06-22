# Population Annealing Monte Carlo (PAMC)

Flexible Population Annealing Monte Carlo (PAMC) framework for general spin models and optimization problems. The code supports 3D Ising and Edwards-Anderson (EA) spin glass models, with adaptive annealing schedules and support for model-specific updates and observables. Designed for research, benchmarking, and pedagogical use, the codebase is structured for extensibility and high-performance simulation workflows.

---

## Project Status

Core functionality is implemented for population annealing on the 3D Ising and EA spin glass models, with adaptive temperature schedules, validated observables, and support for genealogical tracking. The code is structured for modular extension and intended for both research and pedagogical purposes.

A performance-focused refactor is underway to support custom memory pools and compact spin storage (bit arrays), targeting improved efficiency in large-scale spin glass simulations.

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
- [x] Edwards-Anderson spin glass with fully validated benchmark against known ground states

### Infrastructure & Postprocessing

- [x] Unit tests (GoogleTest) and continuous integration
- [ ] Generic observable interface (`measureObservable()`)
- [ ] Data output and I/O framework (`DataWriter` class or equivalent)
- [x] Python postprocessing and analysis scripts for examples and tests
- [x] OpenMP-based parallel update sweeps

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
- `validation/` — Python scripts for validating simulation output and generating analysis plots 
  - `Ising_model/binder_validation.py` — Verify Binder cumulant crossover in 3D Ising model
  - `EA_model/EA_validation.py` — Check PAMC output against known ground state for benchmark disorder realization

---

## Available Models

- 3D Ising model with Metropolis, heat bath, and Wolff updates
- 3D Edwards-Anderson spin glass model with disorder input and benchmark validation

## Planned Features

- Fully connected spin glasses
- Generalized coupling matrices and hybrid optimization support

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
