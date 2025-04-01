# Population Annealing Monte Carlo (PAMC)

Population Annealing Monte Carlo (PAMC) is a flexible and scalable Monte Carlo algorithm for simulating equilibrium and nonequilibrium statistical physics systems in the canonical ensemble. This implementation is designed to be flexible and parallel, supporting OpenMP parallelization, with planned extensions for MPI and GPU acceleration.

## Features
- **Canonical ensemble simulation**  
- **Flexible model support** (starting with the Ising model)  
- **Parallelization options** (initially OpenMP, with future support for MPI and GPUs)  
- **Statistical and systematic error estimators**  
- **Object-oriented C++ implementation** with modular design  

## Installation
### **Prerequisites**
- **C++17 or later**  
- **CMake** (for building the project)  
- **GoogleTest** (for unit testing)  
- **OpenMP** (for parallel execution)  

### **Build Instructions**
Clone the repository and build the project using CMake:
```bash
git clone https://github.com/yourusername/pamc.git
cd pamc
mkdir build && cd build
cmake ..
make -j$(nproc)  # Adjust for your CPU cores