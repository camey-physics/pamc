#ifndef MEMORY_POOL_HPP
#define MEMORY_POOL_HPP

#include <memory>
#include <cassert>
#include <cstddef>

template <typename T>
class MemoryPool {
 public:
  explicit MemoryPool(size_t capacity)
      : capacity_(capacity),
        next_index_(0),
        buffer_(std::make_unique<T[]>(capacity)) {}

  T* allocate(std::size_t count) {
    assert(next_index_ + count <= capacity_);
    T* ptr = &buffer_[next_index_];
    next_index_ += count;
    return ptr;
  }

  void reset() { next_index_ = 0; }

  size_t size() const { return next_index_; } 
  size_t capacity() const { return capacity_; }

 private:
  size_t capacity_;
  size_t next_index_;
  std::unique_ptr<T[]> buffer_;
};

#endif