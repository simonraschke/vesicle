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

#include "monte_carlo.hpp"



void MonteCarlo::setup()
{
    cells.setParameters(getParameters());
    cells.setup();
    cells.deployParticles(*target_range);
    energy_work.reset(new enhance::TriangularMatrix<float>(target_range->size()));
    energy_old.reset(new enhance::TriangularMatrix<float>(target_range->size()));
    sw_position.setup(1000*getParameters().mobile, getParameters().sw_position_min, getParameters().sw_position_max, getParameters().sw_position_target);
    sw_orientation.setup(1000*getParameters().mobile, getParameters().sw_orientation_min, getParameters().sw_orientation_max, getParameters().sw_orientation_target);
}



MonteCarlo::cell_container_type& MonteCarlo::getCells()
{
    return cells;
}



std::size_t MonteCarlo::getCurrentStep() const
{
    return current_step;
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
        const auto num_members_is = cells.membersContained();
        const auto num_members_should = target_range->size();
        if(num_members_is!=num_members_should)
        {
            vesCRITICAL("lost particles while reordering: " << std::boolalpha << (cells.membersContained()==target_range->size()) << " found " << cells.membersContained() << " should be " << target_range->size() )
        }
        assert(num_members_is==num_members_should);
        assert(cells.membersContained()==target_range->size());
        ++current_step;
    }
}



void MonteCarlo::doMCmove(const cell_type& cell)
{
    // vesDEBUG(__PRETTY_FUNCTION__)
    for(Particle& particle : cell)
    {
        // coordinates move
        {
            // const auto stepwidth = getParameters().sw_position;
            const auto stepwidth = sw_position();
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
                // #pragma clang loop vectorize(enable) interleave(enable)
                // for(const Particle& other : region_cell.get())
                // {
                //     if( particle == other ) continue;
                //     (*energy_work)(particle.ID, other.ID) = getInteraction()->potential(particle, other);
                // }
                // #pragma clang loop vectorize(enable) interleave(enable)
                for(std::size_t i = 0; i < region_cell.get().size(); ++i)
                {
                    if( i == particle.ID) continue;
                    (*energy_work)(particle.ID, (region_cell.get().begin()+i)->get().ID) = getInteraction()->potential(particle, (*(region_cell.get().begin()+i)));
                }
                for(const Particle& other : region_cell.get())
                {
                    if( particle == other ) continue;
                    delta_energy += (*energy_work)(particle.ID, other.ID) - (*energy_old)(particle.ID, other.ID);
                }
            }

            // rejection
            if(!acceptance->isValid(delta_energy))
            {
                sw_position.rejected();
                particle.setCoords(particle.coordsOld());
                for(const auto& region_cell : cell.getRegion())
                {
                    // #pragma clang loop vectorize(enable) interleave(enable)
                    for(const Particle& other : region_cell.get())
                    {
                        if( particle == other ) continue;
                        (*energy_work)(particle.ID, other.ID) = (*energy_old)(particle.ID, other.ID);
                    }
                }
            }
            // acctance
            else
                sw_position.accepted();
        }

        // orientation move
        {
            // const float stepwidth = getParameters().sw_orientation;
            const float stepwidth = sw_orientation();
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
                // #pragma clang loop vectorize(enable) interleave(enable)
                for(const Particle& other : region_cell.get())
                {
                    if( particle == other ) continue;
                    (*energy_work)(particle.ID, other.ID) = getInteraction()->potential(particle, other);
                }
                for(const Particle& other : region_cell.get())
                {
                    if( particle == other ) continue;
                    delta_energy += (*energy_work)(particle.ID, other.ID) - (*energy_old)(particle.ID, other.ID);
                }
            }

            // rejection
            if(!acceptance->isValid(delta_energy))
            {
                sw_orientation.rejected();
                particle.setOrientation(particle.orientationOld());
                for(const auto& region_cell : cell.getRegion())
                {
                    // #pragma clang loop vectorize(enable) interleave(enable)
                    for(const Particle& other : region_cell.get())
                    {
                        if( particle == other ) continue;
                        (*energy_work)(particle.ID, other.ID) = (*energy_old)(particle.ID, other.ID);
                    }
                }
            }
            // acctance
            else
                sw_orientation.accepted();
        }
    }
}



float MonteCarlo::getEnergyMatrixSum()
{
    return energy_work->sum();
}