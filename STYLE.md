# C++ Style Guide for PAMC

This style guide defines naming and structural conventions used throughout the Population Annealing Monte Carlo (PAMC) project. The goal is to ensure consistency, readability, and professional-quality code that is both maintainable and extensible.

---

## Naming Conventions

| Element               | Convention       | Example                         |
|-----------------------|------------------|----------------------------------|
| **Classes/Structs**   | `PascalCase`     | `IsingModel`, `SharedModelData` |
| **Methods**           | `camelCase`      | `initializeState()`, `updateSweep()` |
| **Member Variables**  | `snake_case_`    | `spin_vector_`, `num_spins_`    |
| **Local Variables**   | `snake_case`     | `system_size`, `replica_index`  |
| **Constants**         | `ALL_CAPS` (with `constexpr`) | `constexpr int MAX_STEPS = 100;` |
| **Enums (scoped)**    | `PascalCase` with `SCREAMING_CASE` values | `enum class UpdateMethod { SWEEP, CLUSTER };` |
| **Template Parameters** | `PascalCase`   | `template <typename ModelT>`    |

---

## File Naming

| Type                  | Convention       | Example              |
|-----------------------|------------------|----------------------|
| Header Files          | `PascalCase.hpp` | `IsingModel.hpp`     |
| Implementation Files  | `PascalCase.cpp` | `IsingModel.cpp`     |
| Shared Templates      | `PascalCase.hpp` | `SharedModelData.hpp`|

---

## File Organization

- **`include/`**: Public headers, abstract interfaces (e.g., `Model.hpp`, `Population.hpp`)
- **`include/models/`**: Model headers (e.g., `IsingModel.hpp`)
- **`src/`**: Core algorithm components (e.g., `Population.cpp`)
- **`src/models/`**: Model implementations (e.g., `IsingModel.cpp`)
- **`tests/`**: Unit tests (GoogleTest)

---

## Code Structure Guidelines

- Always prefer `const` where possible
- Pass large objects by reference or pointer
- Use `std::span<T>` when passing views into shared or contiguous data
- Avoid raw pointers unless for non-owning access to external memory (e.g., pooled data)

---

## Comments

- Use `//` for inline or brief comments
- Use `///` or `/** */` for Doxygen-style documentation (optional)
- Explain reasoning, not just what code is doing

---

## Future Tooling

To support automated formatting, a `.clang-format` file may be added to enforce these conventions.

---

This style guide applies to all components in the PAMC codebase. For consistency, contributors and future refactors should adhere to these guidelines unless there's a strong reason to deviate.