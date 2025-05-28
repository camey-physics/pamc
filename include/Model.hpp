#ifndef MODEL_HPP
#define MODEL_HPP

class Model {
 public:
  virtual ~Model() = default;
  virtual void initializeState() = 0;
  virtual void copyStateFrom(const Model& other) = 0;

  enum class UpdateMethod {};
  virtual double calcEnergy() const = 0;
  virtual void updateSweep(int num_sweeps) = 0;
};

#endif  // MODEL_HPP