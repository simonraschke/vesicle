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

#include "data_collector.hpp"



DataCollector::DataCollector()
    : working_dir(boost::filesystem::current_path())
{
}



void DataCollector::setup()
{
    if(boost::filesystem::exists(getParameters().analysis_path) && !getParameters().analysis_overwrite)
    {
        FILE = std::make_unique<HighFive::File>
        (
            boost::filesystem::system_complete(working_dir/getParameters().analysis_path).string(), 
            HighFive::File::ReadWrite
        );
        FILE->createGroup("/cluster_histograms");
    }
    else
    {
        FILE = std::make_unique<HighFive::File>
        (
            boost::filesystem::system_complete(working_dir/getParameters().analysis_path).string(), 
            HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate
        );
        FILE->createGroup("/cluster_histograms");
    }
    for(const auto& [key,val] : systemAttributes())
    {
        HighFive::Attribute attribute = FILE->getGroup("/cluster_histograms").createAttribute<std::string>(key, HighFive::DataSpace::From(val));
        attribute.write(val);
    }

    if(enhance::splitAtDelimiter(getParameters().analysis_input.string(),".").back() == "gro")
        reader = std::make_unique<TrajectoryReaderGro>();
    reader->setParameters(getParameters());
    reader->setPath(getParameters().analysis_input);

    clusters.setParameters(getParameters());
}



void DataCollector::collect()
{
    while(!reader->isEOF() && Controller::SIGNAL == 0)
    {
        reader->readNextFrame(getParameters().analysis_frames);
        try
        {
            auto snapshot = reader->getFrame(-1);
            translator(snapshot);
            
            timepoints.emplace_back(translator.time_elapsed);

            if(getParameters().analysis_full || getParameters().analysis_epot)
                try_potential_energy();


            if(getParameters().analysis_full || getParameters().analysis_cluster)
                try_cluster();
        }
        catch (const std::exception& e)
        {
            vesWARNING(e.what())
        }
        vesLOG( "EOF: " << std::boolalpha << reader->isEOF() )
    }
}



void DataCollector::write()
{
    tbb::parallel_invoke
    (
        // potential energy
        [&]
        {
            assert(timepoints.size() == potential_energies.size());
            boost::multi_array<float,2> potential_energies_array(boost::extents[potential_energies.size()][2]);
            for(std::size_t i = 0; i < potential_energies.size(); ++i)
            {
                potential_energies_array[i][0] = timepoints[i];
                potential_energies_array[i][1] = potential_energies[i];
            }

            HighFive::DataSet dataset = FILE->createDataSet<float>("/potential_energies", HighFive::DataSpace::From(potential_energies_array));
            for(const auto& [key,val] : systemAttributes())
            {
                HighFive::Attribute attribute = dataset.createAttribute<std::string>(key, HighFive::DataSpace::From(val));
                attribute.write(val);
            }
            dataset.write(potential_energies_array);
        },

        // clusters volumes
        [&]
        {
            assert(timepoints.size() == cluster_volumes.size());
            boost::multi_array<float,2> cluster_volumes_array(boost::extents[cluster_volumes.size()][2]);
            for(std::size_t i = 0; i < cluster_volumes.size(); ++i)
            {
                cluster_volumes_array[i][0] = timepoints[i];
                cluster_volumes_array[i][1] = cluster_volumes[i];
            }

            HighFive::DataSet dataset = FILE->createDataSet<float>("/cluster_volumes", HighFive::DataSpace::From(cluster_volumes_array));
            for(const auto& [key,val] : systemAttributes())
            {
                HighFive::Attribute attribute = dataset.createAttribute<std::string>(key, HighFive::DataSpace::From(val));
                attribute.write(val);
            }
            dataset.write(cluster_volumes_array);
        },

        // clusters surface areas
        [&]
        {
            assert(timepoints.size() == cluster_surface_areas.size());
            boost::multi_array<float,2> cluster_surface_areas_array(boost::extents[cluster_surface_areas.size()][2]);
            for(std::size_t i = 0; i < cluster_surface_areas.size(); ++i)
            {
                cluster_surface_areas_array[i][0] = timepoints[i];
                cluster_surface_areas_array[i][1] = cluster_surface_areas[i];
            }

            HighFive::DataSet dataset = FILE->createDataSet<float>("/cluster_surface_areas", HighFive::DataSpace::From(cluster_surface_areas_array));
            for(const auto& [key,val] : systemAttributes())
            {
                HighFive::Attribute attribute = dataset.createAttribute<std::string>(key, HighFive::DataSpace::From(val));
                attribute.write(val);
            }
            dataset.write(cluster_surface_areas_array);
        }
    );
}



void DataCollector::try_potential_energy()
{
    PotentialEnergyParser parser;
    parser.setParameters(getParameters());
    parser.setTarget(&translator.particles);

    if(getParameters().interaction == "alj")
        parser.setInteraction<AngularLennardJones>();
    else if(getParameters().interaction == "lj")
        parser.setInteraction<LennardJones>();

    parser.parse();

    potential_energies.emplace_back(parser.result);
}



void DataCollector::try_cluster()
{
    {
        clusters.setTarget(translator.particles.begin(),translator.particles.end());
        if(getParameters().analysis_cluster_algorithm == "DBSCAN")
            clusters.DBSCANrecursive(getParameters().analysis_cluster_minimum_size, getParameters().analysis_cluster_distance_threshold);
        else
            vesCRITICAL("cluster algorithm " << getParameters().analysis_cluster_algorithm << " unknown")
    }

    {
        boost::multi_array<long int,2> cluster_histogram(boost::extents[clusters.numClusters()][clusters.maxClusterSize()+1]);

        for(std::size_t i = 0; i != clusters.numClusters(); ++i) 
        for(std::size_t j = 0; j != clusters.maxClusterSize()+1; ++j)
            cluster_histogram[i][j] = -1;

        for(std::size_t i = 0; i < clusters.numClusters(); ++i)
        {
            cluster_histogram[i][0] = (clusters.begin()+i)->size();
            for(std::size_t j = 0; j < (clusters.begin()+i)->size(); ++j)
            {
                cluster_histogram[i][j+1] = (*((clusters.begin()+i)->begin()+j))->ID;
            }
        }

        HighFive::DataSet dataset = FILE->createDataSet<long int>("/cluster_histograms/time"+boost::lexical_cast<std::string>(timepoints.back()), HighFive::DataSpace::From(cluster_histogram));
        dataset.write(cluster_histogram);
    }

    {
        tbb::concurrent_vector<ClusterStructureParser> structureParsers;
        tbb::parallel_for_each(clusters.begin(), clusters.end(), [&](const auto& cluster)
        {
            if(cluster.size() >= getParameters().analysis_cluster_significant_size)
            {
                auto it = structureParsers.emplace_back(cluster);
                it->setParameters(getParameters());
                it->parse();
            }
        });

        tbb::parallel_invoke
        (
            [&]
            { 
                cluster_volumes.emplace_back(PARALLEL_REDUCE(float, structureParsers, [&](float i, const auto& parser)
                {
                    return i + parser.getVolume();
                }));
            },

            [&]
            { 
                cluster_surface_areas.emplace_back(PARALLEL_REDUCE(float, structureParsers, [&](float i, const auto& parser)
                {
                    return i + parser.getSurfaceArea();
                }));
            },
            
            [&]
            { 
                if(!structureParsers.empty())
                std::max_element(std::begin(structureParsers), std::end(structureParsers), [](auto& a, auto& b)
                {
                    return a.result < b.result;
                })->printXML("largest.vtu");
            }
        );
    }
}