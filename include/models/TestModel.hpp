#ifndef TESTMODEL_HPP
#define TESTMODEL_HPP


#include "models/Model.hpp"

class TestModel : public Model {
    public:
        TestModel() = default;
        explicit TestModel(double state);
        // ~TestModel();
        void initializeState() override;
        void initializeState(double state);
        void copyStateFrom(const Model& other) override;
        enum class UpdateMethod { };
        double calcEnergy() const override;
        void updateSweep(int numSweeps) override;
        // I am not sure about creating new clones as smart pointers. I am not sure that this will be the most
        // efficient way to deal with replicas and probably prefer to do this within a separate population class
        // virtual std::unique_ptr<Model> clone() const = 0; 
    private:
        double state_;
    };

    #endif // TESTMODEL_HPP