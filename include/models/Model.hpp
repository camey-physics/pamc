#ifndef MODEL_HPP
#define MODEL_HPP

class Model {
public:
    virtual ~Model() = default;
    virtual void initializeState() = 0;
    
    virtual void copyStateFrom(const Model& other) = 0;
    enum class UpdateMethod { };
    virtual double calcEnergy() const = 0;
    virtual void updateSweep(int numSweeps) = 0;
    // I am not sure about creating new clones as smart pointers. I am not sure that this will be the most
    // efficient way to deal with replicas and probably prefer to do this within a separate population class
    // virtual std::unique_ptr<Model> clone() const = 0; 

};

#endif // MODEL_HPP