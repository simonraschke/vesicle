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

#include "cluster_parser.hpp"



// template<PERIODIC P>
// void ClusterParser<P>::setTarget(PARTICLERANGE* range)
// {
//     vesDEBUG(__PRETTY_FUNCTION__)
//     target_range = enhance::make_observer<PARTICLERANGE>(range);
//     particles.clear();
//     clusters.clear();
//     tbb::parallel_for_each(std::begin(*range), std::end(*range),[&](const PARTICLERANGE::value_type& particle)
//     {
//         particles.emplace_back(particle.get());
//     });


//     vesLOG("DBSCAN analysing " << particles.size() << " particles")
// }



