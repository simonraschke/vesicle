#include "particle_distributor.hpp"



void RandomDistributor::operator()(PARTICLERANGE* range)
{
    assert(range);
    tbb::parallel_for_each(range->begin(), range->end(), [&](auto& p) 
    {
        assert(p);
        p->setCoords(randomCoords());
    });

    assert(range);
    int i = 0;
    for(auto& p : *range)
    {
        std::cout << i++ << '\n';
        while(conflicting_placement(range,p))
        {
            assert(p);
            p->setCoords(randomCoords());
        }
    }
}



RandomDistributor::cartesian RandomDistributor::randomCoords() const
{
    return cartesian
    (
        enhance::random<cartesian::Scalar>(0.f,getLengthX()),
        enhance::random<cartesian::Scalar>(0.f,getLengthY()),
        enhance::random<cartesian::Scalar>(0.f,getLengthZ())
    );
}



bool RandomDistributor::conflicting_placement(PARTICLERANGE* range, PARTICLERANGE::value_type& p1)
{
    std::atomic<bool> conflict {false};
    assert(range);
    tbb::parallel_for_each(range->begin(), range->end(), [&](auto& p2) 
    {
        assert(p1 && p2);
        if(p1==p2) return;
        if(squared_distance(*p1,*p2) <= 1.f) conflict.store(true);
    });
    return conflict.load();
}