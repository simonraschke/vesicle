/*  
*   Copyright 2017 Simon Raschke
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

#include "monte_carlo.hpp"



void MonteCarlo::setup()
{
    cells.setParameters(getParameters());
    cells.setup();
    cells.deployParticles(*target_range);
    energy_work.reset(new enhance::TriangularMatrix<float>(target_range->size()));
    energy_old.reset(new enhance::TriangularMatrix<float>(target_range->size()));
}



void MonteCarlo::step(const unsigned long& steps)
{
    vesDEBUG(__PRETTY_FUNCTION__)
    for(unsigned long step = 0; step < steps; ++step)
    {
        assert(energy_old);
        assert(energy_work);
        std::swap(energy_old, energy_work);
        assert(target_range);
        tbb::parallel_for_each(target_range->begin(), target_range->end(), [&](auto& target) 
        {
            assert(target);
            target->save();
        });

        CellBasedAlgorithm::preparation(cells);

        CellBasedAlgorithm::step(cells, [&](const Cell<Particle>& cell){ doMCmove(cell);});

        assert(cells.membersContained()==target_range->size());
        CellBasedAlgorithm::reorder(cells);
        vesDEBUG("lost particles while reordering: " << std::boolalpha << (cells.membersContained()==target_range->size()) << " found " << cells.membersContained() << " should be " << target_range->size() )
        assert(cells.membersContained()==target_range->size());
    }
}




// float MonteCarlo::potentialEnergy(const std::unique_ptr<Particle>& p1) const
// {
//     std::atomic<float> energy_sum {0.f};
//     tbb::parallel_for_each(target_range->begin(), target_range->end(), [&](auto& target) 
//     {
//         assert(p1);
//         assert(target);
//         if(p1==target) return;
//         assert(getInteraction());
//         const float energy = getInteraction()->potential(p1,target);

//         auto current = energy_sum.load();
//         while (!energy_sum.compare_exchange_weak(current, current + energy))
//             current = energy_sum.load();
//     });
//     return !std::isnan(energy_sum.load()) ? energy_sum.load() : throw std::runtime_error("potential Energy is NAN");
// }



// float MonteCarlo::potentialEnergyInRegion(const cell_type& cell, const Particle& particle) const
// {
//     // vesDEBUG(__PRETTY_FUNCTION__)
//     float energy_sum = 0.f;
//     for(const cell_type& other_cell: cell.getRegion())
//     for(const Particle& other: other_cell)
//     // std::for_each(cell.getRegion().begin(), cell.getRegion().end(), [&](const Particle& other) 
//     {
//         if(std::addressof(particle)==std::addressof(other)) continue;
//         assert(getInteraction());
//         energy_sum += getInteraction()->potential(particle,other);
//     }
//     return !std::isnan(energy_sum) ? energy_sum : throw std::runtime_error("potential Energy is NAN");
// }



void MonteCarlo::doMCmove(const cell_type& cell)
{
    // vesDEBUG(__PRETTY_FUNCTION__)
    for(Particle& particle : cell)
    {
        // coordinates move
        {
            const auto stepwidth = getParameters().stepwidth_coordinates;
            const auto translation = Particle::cartesian
            (
                enhance::random<float>(-stepwidth,stepwidth),
                enhance::random<float>(-stepwidth,stepwidth),
                enhance::random<float>(-stepwidth,stepwidth)
            );

            particle.setCoords(particle.coords()+translation);
            
            float delta_energy = 0.f;
            for(const auto& region_cell : cell.getRegion())
            {
                for(const Particle& other : region_cell.get())
                {
                    if( particle == other ) continue;
                    (*energy_work)(particle.ID, other.ID) = getInteraction()->potential(particle, other);
                    delta_energy += (*energy_work)(particle.ID, other.ID) - (*energy_old)(particle.ID, other.ID);
                }
            }

            // rejection
            if(!acceptance->isValid(delta_energy))
            {
                particle.setCoords(particle.coordsOld());
                for(const auto& region_cell : cell.getRegion())
                {
                    for(const Particle& other : region_cell.get())
                    {
                        if( particle == other ) continue;
                        (*energy_work)(particle.ID, other.ID) = (*energy_old)(particle.ID, other.ID);
                    }
                }
            }
        }

        // orientation move
        {
            const float stepwidth = getParameters().stepwidth_orientation;
            const auto orientation = Particle::cartesian
            (
                enhance::random<float>(-1.f,1.f),
                enhance::random<float>(-1.f,1.f),
                enhance::random<float>(-1.f,1.f)
            );
            const Eigen::AngleAxisf rotation (stepwidth, orientation);
            particle.setOrientation(rotation * particle.orientation());

            float delta_energy = 0.f;
            for(const auto& region_cell : cell.getRegion())
            {
                for(const Particle& other : region_cell.get())
                {
                    if( particle == other ) continue;
                    (*energy_work)(particle.ID, other.ID) = getInteraction()->potential(particle, other);
                    delta_energy += (*energy_work)(particle.ID, other.ID) - (*energy_old)(particle.ID, other.ID);
                }
            }
            // const float energy_before = potentialEnergyInRegion(cell,particle);
            // const float energy_after = potentialEnergyInRegion(cell,particle);

            // rejection
            if(!acceptance->isValid(delta_energy))
            {
                particle.setOrientation(particle.orientationOld());
                for(const auto& region_cell : cell.getRegion())
                {
                    for(const Particle& other : region_cell.get())
                    {
                        if( particle == other ) continue;
                        (*energy_work)(particle.ID, other.ID) = (*energy_old)(particle.ID, other.ID);
                    }
                }
            }
        }
    }
}
