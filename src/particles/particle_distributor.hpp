#pragma once

#include "definitions.hpp"
#include "particle.hpp"
#include "vesicleIO/parameters.hpp"
#include "enhance/random.hpp"
#include "systems/box.hpp"



struct Distributor
    : public Box<PERIODIC::ON>
    , virtual public ParameterDependentComponent
{
    virtual void operator()(PARTICLERANGE*) = 0;
    
    virtual ~Distributor() = default;

protected:
    Distributor() = default;
};



struct RandomDistributor
    : public Distributor
{
    typedef PARTICLERANGE::value_type::element_type::cartesian cartesian;

    virtual void operator()(PARTICLERANGE* range) override
    {
        for(const auto& p : *range)
        {
            {
                cartesian vec;
                vec(0) = enhance::random<cartesian::Scalar>(0.f,getLengthX());
                vec(1) = enhance::random<cartesian::Scalar>(0.f,getLengthY());
                vec(2) = enhance::random<cartesian::Scalar>(0.f,getLengthZ());
                assert(p);
                p->updateCoords(vec);
            }
            {
                cartesian vec = p->coords();
                vec(0) += enhance::random<cartesian::Scalar>(-getParameters().dt,getParameters().dt);
                vec(1) += enhance::random<cartesian::Scalar>(-getParameters().dt,getParameters().dt);
                vec(2) += enhance::random<cartesian::Scalar>(-getParameters().dt,getParameters().dt);
                assert(p);
                p->updateCoords(vec);
            }
        }
    }
};