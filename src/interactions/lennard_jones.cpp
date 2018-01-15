/*  
*   Copyright 2017-2018 Simon Raschke
*
*   Licensed under the Apache License, Version 2.0 (the "License");
*   you may not use this file except in compliance with the License.
*   You may obtain a copy of the License at
*
*       http://www.apache.org/licenses/LICENSE-2.0
*
*   Unless required by applicable law or agreed to in writing, software
*   distributed under the License is distributed on an "AS IS" BASIS,
*   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*   See the License for the specific language governing permissions and
*   limitations under the License.
*/

#include "lennard_jones.hpp"



float LennardJones::potential(const Particle& p1, const Particle& p2) const 
{
    const float r2 = 1.f/squared_distance(p1,p2);
    const float r6 = r2*r2*r2;
    return 4.f*(r6*r6-r6);
}



LennardJones::cartesian LennardJones::force(const Particle& p1, const Particle& p2) const 
{
    const float r2 = 1.f/squared_distance(p1.coordsOld(),p2.coordsOld());
    const float r6 = r2*r2*r2;
    const float value = -24.f*r2*r6*(r6*2-1.f);
#ifndef NDEBUG
    if(!std::isfinite(squared_distance(p1,p2)) || !std::isfinite(r2) || !std::isfinite(value) || distanceVector(p1,p2).hasNaN()) 
    {
        vesWARNING("distanceVector " << distanceVector(p1,p2).format(ROWFORMAT))
        vesWARNING("squared_distance(p1,p2) " << squared_distance(p1,p2))
        vesWARNING("r2 " << r2)
        vesWARNING("r6 " << r6)
        vesWARNING("value " << value)
    }
#endif

    return distanceVector(p1.coordsOld(),p2.coordsOld())*value;
}



bool LennardJones::isAnisotropic() const
{
    return false;
}