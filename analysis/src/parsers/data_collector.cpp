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
        vesLOG("HDF5: create group /cluster_frame_guided")
        FILE->createGroup("/cluster_frame_guided");
        vesLOG("HDF5: create group /cluster_self_assembled")
        FILE->createGroup("/cluster_self_assembled");
    }
    else
    {
        FILE = std::make_unique<HighFive::File>
        (
            boost::filesystem::system_complete(working_dir/getParameters().analysis_path).string(), 
            HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate
        );
        vesLOG("HDF5: create group /cluster_frame_guided")
        FILE->createGroup("/cluster_frame_guided");
        vesLOG("HDF5: create group /cluster_self_assembled")
        FILE->createGroup("/cluster_self_assembled");
    }
    for(const auto& [key,val] : systemAttributes())
    {
        {
            HighFive::Attribute attribute = FILE->getGroup("/cluster_frame_guided").createAttribute<std::string>(key, HighFive::DataSpace::From(val));
            vesLOG("HDF5: write attributes to group /cluster_frame_guided")
            attribute.write(val);
        }
        {
            HighFive::Attribute attribute = FILE->getGroup("/cluster_self_assembled").createAttribute<std::string>(key, HighFive::DataSpace::From(val));
            vesLOG("HDF5: write attributes to group /cluster_self_assembled")
            attribute.write(val);
        }
    }
    {
        std::string size("size");
        std::string order("order");
        std::string volume("volume");
        std::string surfacearea("surfacearea");
        HighFive::Attribute a = FILE->getGroup("/cluster_frame_guided").createAttribute<std::string>("col0", HighFive::DataSpace::From(size));
        HighFive::Attribute b = FILE->getGroup("/cluster_frame_guided").createAttribute<std::string>("col1", HighFive::DataSpace::From(order));
        HighFive::Attribute c = FILE->getGroup("/cluster_frame_guided").createAttribute<std::string>("col2", HighFive::DataSpace::From(volume));
        HighFive::Attribute d = FILE->getGroup("/cluster_frame_guided").createAttribute<std::string>("col3", HighFive::DataSpace::From(surfacearea));
        HighFive::Attribute e = FILE->getGroup("/cluster_self_assembled").createAttribute<std::string>("col0", HighFive::DataSpace::From(size));
        HighFive::Attribute f = FILE->getGroup("/cluster_self_assembled").createAttribute<std::string>("col1", HighFive::DataSpace::From(order));
        HighFive::Attribute g = FILE->getGroup("/cluster_self_assembled").createAttribute<std::string>("col2", HighFive::DataSpace::From(volume));
        HighFive::Attribute h = FILE->getGroup("/cluster_self_assembled").createAttribute<std::string>("col3", HighFive::DataSpace::From(surfacearea));
        a.write(size);
        b.write(order);
        c.write(volume);
        d.write(surfacearea);
        e.write(size);
        f.write(order);
        g.write(volume);
        h.write(surfacearea);
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
            vesCRITICAL(e.what())
        }
        vesLOG( "EOF: " << std::boolalpha << reader->isEOF() )
    }
}



void DataCollector::write()
{
    // potential energy
    {
        assert(timepoints.size() == potential_energies.size());
        boost::multi_array<float,2> potential_energies_array(boost::extents[potential_energies.size()][2]);
        for(std::size_t i = 0; i < potential_energies.size(); ++i)
        {
            potential_energies_array[i][0] = timepoints[i];
            potential_energies_array[i][1] = potential_energies[i];
        }

        vesLOG("HDF5: create dataset /potential_energies")
        HighFive::DataSet dataset = FILE->createDataSet<float>("/potential_energies", HighFive::DataSpace::From(potential_energies_array));
        for(const auto& [key,val] : systemAttributes())
        {
            HighFive::Attribute attribute = dataset.createAttribute<std::string>(key, HighFive::DataSpace::From(val));
            attribute.write(val);
        }
        vesLOG("HDF5: write dataset /potential_energies with size " << potential_energies.size() << "x2")
        dataset.write(potential_energies_array);
    }
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
    // perform cluster algorithm
    {
        clusters.setTarget(translator.particles.begin(),translator.particles.end());
        clusters.scan();
    }

    {
        // structureParsers contains all cluster of significant size
        // every parser references its origin cluster
        // tbb::concurrent_vector<ClusterStructureParser> structureParsers;
        // for powercrust ONLY single threaded!
        tbb::parallel_for_each(clusters.begin(), clusters.end(), [&](auto& cluster)
        {
            cluster.rearrangeSubclusters();
        });

        vesLOG(__PRETTY_FUNCTION__)

        // write histograms every step
        write_fg_histogram(); 
        write_sa_histogram(); 

        // write a structure file of the largest cluster
        { 
            if(!clusters.empty())
            {
                std::max_element(std::begin(clusters), std::end(clusters), [](auto& a, auto& b)
                {
                    return a.getStructure().getVolume() < b.getStructure().getVolume();
                })->getStructure().printXML("largest.vtu");
            }
        }
    }
}


void DataCollector::write_fg_histogram()
{
    typedef float hist_t;
    const auto largestFGCluster = clusters.getLargest([](const auto& cluster)
    { 
        return ! (cluster.template containsParticleType<PARTICLETYPE::FRAME>()); 
    });

    const std::size_t numFrameGuidedClusters = std::count_if(std::begin(clusters), std::end(clusters), [](const auto& cluster)
    { 
        return cluster.template containsParticleType<PARTICLETYPE::FRAME>();
    });

    if(!largestFGCluster || numFrameGuidedClusters == 0)
        return;

    const auto maxClusterFrameGuided = largestFGCluster->size();

    boost::multi_array<hist_t,2> cluster_histogram(boost::extents[numFrameGuidedClusters][maxClusterFrameGuided+4]);

    // initialize array with -1
    for(std::size_t i = 0; i != numFrameGuidedClusters; ++i) 
        for(std::size_t j = 0; j != maxClusterFrameGuided+4; ++j)
            cluster_histogram[i][j] = static_cast<hist_t>(-1);

    // matrix with 
    // column 0 size of cluster
    // further columns with particle IDs in this cluster
    for(std::size_t i = 0; i < numFrameGuidedClusters; ++i)
    {   
        const auto& cluster = *(clusters.begin()+i);
        if( !cluster.containsParticleType<PARTICLETYPE::FRAME>() )
            continue;

        // size
        cluster_histogram[i][0] = static_cast<hist_t>(cluster.size());
        cluster_histogram[i][1] = static_cast<hist_t>(cluster.getOrder());
        cluster_histogram[i][2] = static_cast<hist_t>(cluster.getStructure().getVolume());
        cluster_histogram[i][3] = static_cast<hist_t>(cluster.getStructure().getSurfaceArea());

        // all particles of cluster
        for(std::size_t j = 0; j < cluster.size(); ++j)
        {
            const auto& particle = *(cluster.begin()+j);

            // ID of particle in cluster
            cluster_histogram[i][j+4] = static_cast<hist_t>(particle.ID);
        }
    }

    vesLOG("HDF5: create dataset /cluster_frame_guided/time"+boost::lexical_cast<std::string>(timepoints.back()))
    HighFive::DataSet dataset = FILE->createDataSet<hist_t>("/cluster_frame_guided/time"+boost::lexical_cast<std::string>(timepoints.back()), HighFive::DataSpace::From(cluster_histogram));
    vesLOG("HDF5: write dataset  /cluster_frame_guided/time"+boost::lexical_cast<std::string>(timepoints.back())+" with size " << numFrameGuidedClusters << "x" << maxClusterFrameGuided+1)
    dataset.write(cluster_histogram);
}



void DataCollector::write_sa_histogram()
{
    typedef float hist_t;
    const std::size_t numSelfAssembledClusters = std::count_if(std::begin(clusters), std::end(clusters),[](const auto& cluster){ return ! (cluster.template containsParticleType<PARTICLETYPE::FRAME>()); });
    const std::size_t maxClusterSelfAssembled = clusters.getLargest().size();
    
    boost::multi_array<hist_t,2> cluster_histogram(boost::extents[numSelfAssembledClusters][maxClusterSelfAssembled+4]);

    if(numSelfAssembledClusters == 0)
        return;
    else if(maxClusterSelfAssembled == 0)
        throw std::logic_error("the max self-assembled cluster is of size 0");

    // initialize array with -1
    for(std::size_t i = 0; i != numSelfAssembledClusters; ++i) 
        for(std::size_t j = 0; j != maxClusterSelfAssembled+4; ++j)
            cluster_histogram[i][j] = static_cast<hist_t>(-1);

    // matrix with 
    // column 0 size of cluster
    // further columns with particle IDs in this cluster
    for(std::size_t i = 0; i < numSelfAssembledClusters; ++i)
    {   
        const auto& cluster = *(clusters.begin()+i);
        if( cluster.containsParticleType<PARTICLETYPE::FRAME>() )
            continue;

        // size
        cluster_histogram[i][0] = static_cast<hist_t>(cluster.size());
        cluster_histogram[i][1] = static_cast<hist_t>(cluster.getOrder());
        cluster_histogram[i][2] = static_cast<hist_t>(cluster.getStructure().getVolume());
        cluster_histogram[i][3] = static_cast<hist_t>(cluster.getStructure().getSurfaceArea());

        // all particles of cluster
        for(std::size_t j = 0; j < cluster.size(); ++j)
        {
            const auto& particle = *(cluster.begin()+j);

            // ID of particle in cluster
            cluster_histogram[i][j+4] = static_cast<hist_t>(particle.ID);
        }
    }

    vesLOG("HDF5: create dataset /cluster_self_assembled/time"+boost::lexical_cast<std::string>(timepoints.back()))
    HighFive::DataSet dataset = FILE->createDataSet<hist_t>("/cluster_self_assembled/time"+boost::lexical_cast<std::string>(timepoints.back()), HighFive::DataSpace::From(cluster_histogram));
    vesLOG("HDF5: write dataset  /cluster_self_assembled/time"+boost::lexical_cast<std::string>(timepoints.back())+" with size " << numSelfAssembledClusters << "x" << maxClusterSelfAssembled+1);
    dataset.write(cluster_histogram);
}
