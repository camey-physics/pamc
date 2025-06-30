#ifndef MEMORY_BLOCK_HPP
#define MEMORY_BLOCK_HPP

// MemoryBlock is used by models to specify memory requirements for external MemoryPool allocation.

#include <cstddef>
#include <typeindex>
#include <typeinfo>

struct MemoryBlock {
  std::type_index type;
  std::size_t count;

  template <typename T>
  static MemoryBlock forType(std::size_t count) {
    return { std::type_index(typeid(T)), count };
  }
};

#endif