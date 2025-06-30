
#include <typeindex>

struct MemoryBlock {
  std::type_index type;
  std::size_t count;
};