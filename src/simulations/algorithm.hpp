#pragma once

#include "definitions.hpp"
#include "enhance/observer_ptr.hpp"
#include "vesicleIO/parameters.hpp"
#include "particles/particle.hpp"
#include "interactions/interaction.hpp"


class Algorithm
    : public ParameterDependentComponent
{
public:

    // set Parameters
    void setTarget(PARTICLERANGE*);

    template<typename I>
    void setInteraction();
    Interaction& getInteraction() const;

    // execute
    virtual void step(const unsigned long& = 1) = 0;

    // necessary
    virtual ~Algorithm() = default;

protected:
    Algorithm() = default;

    virtual void updateCoords() = 0;
    virtual void updateOrientations() = 0;
    virtual void updateForces() = 0;
    virtual void updateVelocities() = 0;

    std::unique_ptr<Interaction> interaction {nullptr};
    enhance::observer_ptr<PARTICLERANGE> target_range {nullptr};

private:  

};



template<typename I>
void Algorithm::setInteraction()
{
    interaction = std::make_unique<I>();
    interaction->setParameters(getParameters());
    interaction->setup();
}