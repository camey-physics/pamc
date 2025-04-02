#include "models/Model.hpp"

class TestModel : public Model {
    public:
        explicit TestModel(double state);
        ~TestModel();
        void initializeState() override;
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