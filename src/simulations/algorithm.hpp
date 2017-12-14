#pragma once

#include "definitions.hpp"
#include "enhance/observer_ptr.hpp"
#include "vesicleIO/parameters.hpp"
#include "particles/particle.hpp"
#include "interactions/interaction.hpp"
#include "acceptance_adapters/acceptance_adapter.hpp"
#include <tbb/tbb.h>
#include <eigen3/Eigen/Sparse>


class Algorithm
    : public ParameterDependentComponent
{
public:

    // set Parameters
    void setTarget(PARTICLERANGE*);

    template<typename I>
    void setInteraction();
    const std::unique_ptr<Interaction>& getInteraction() const;

    // this should be nullptr if not MonteCarlo
    template<typename A>
    void setAcceptance();
    const std::unique_ptr<AcceptanceAdapter>& getAcceptance() const;

    // execute
    virtual void step(const unsigned long& = 1) = 0;

    // necessary
    virtual ~Algorithm() = default;

protected:
    Algorithm() = default;

    virtual void updateCoords() = 0;
    // virtual void updateOrientations() = 0;
    virtual void updateForces() = 0;
    virtual void updateVelocities() = 0;

    std::unique_ptr<Interaction> interaction {nullptr};
    std::unique_ptr<AcceptanceAdapter> acceptance {nullptr};
    enhance::observer_ptr<PARTICLERANGE> target_range {nullptr};

    // std::unique_ptr<Eigen::SparseMatrix<float>> energy_matrix_old  {std::make_unique<Eigen::SparseMatrix<float>>()};
    // std::unique_ptr<Eigen::SparseMatrix<float>> energy_matrix_work {std::make_unique<Eigen::SparseMatrix<float>>()};

private:  

};



template<typename I>
void Algorithm::setInteraction()
{
    vesDEBUG(__PRETTY_FUNCTION__)
    interaction = std::make_unique<I>();
    interaction->setParameters(getParameters());
    interaction->setup();
}



template<typename A>
void Algorithm::setAcceptance()
{
    vesDEBUG(__PRETTY_FUNCTION__)
    acceptance = std::make_unique<A>();
    acceptance->setParameters(getParameters());
    // acceptance->setup();
}