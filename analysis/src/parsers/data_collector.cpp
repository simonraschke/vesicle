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

    if(enhance::splitAtDelimiter(getParameters().analysis_input.string(),".").back() == "gro")
        reader = std::make_unique<TrajectoryReaderGro>();
    reader->setParameters(getParameters());
    reader->setPath(getParameters().analysis_input);

    clusters.setParameters(getParameters());
}



void DataCollector::collect()
{
    // // we create a new hdf5 file
    // const unsigned long x = 100000000;
    // const unsigned long y = 2;

    // boost::multi_array<float, 2> my_array(boost::extents[x][y]);
    // for (size_t i = 0; i < x; ++i) {
    //     for (size_t j = 0; j < y; ++j) {
    //         my_array[i][j] = 1.337 + i + j;
    //     }
    // }

    // HighFive::File file(boost::filesystem::system_complete("data.h5").string(), HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate);
    // HighFive::DataSet dataset = file.createDataSet<float>("/boost_multi_array", HighFive::DataSpace::From(my_array));

    // // lets write our vector of int to the HDF5 dataset
    // dataset.write(my_array);

    // read back
    // boost::multi_array<int, 2> result;
    // dataset.read(result);

    // for(auto i : result)
    // for(auto j : i)
    // vesLOG(j)

    while(!reader->isEOF() && Controller::SIGNAL == 0)
    {
        reader->readNextFrame(std::regex("^[0-9]*$"));
        try
        {
            auto frame = reader->getFrame(-1);
            translator(frame);
            
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
    // potential energy
    {
        assert(timepoints.size() == potential_energies.size());
        boost::multi_array<float,2> potential_energies_array(boost::extents[potential_energies.size()][2]);
        for(std::size_t i = 0; i < potential_energies.size(); ++i)
        {
            potential_energies_array[i][0] = timepoints[i];
            potential_energies_array[i][1] = potential_energies[i];
        }

        HighFive::DataSet dataset = FILE->createDataSet<float>("/potential_energies", HighFive::DataSpace::From(potential_energies_array));
        dataset.write(potential_energies_array);
    }

    // clusters
    // {
    //     assert(timepoints.size() == potential_energies.size());
    //     boost::multi_array<float,2> potential_energies_array(boost::extents[translator.particles.size()][2]);
    //     for(std::size_t i = 0; i < potential_energies.size(); ++i)
    //     {
    //         potential_energies_array[i][0] = timepoints[i];
    //         potential_energies_array[i][1] = potential_energies[i];
    //     }

    //     HighFive::DataSet dataset = FILE->createDataSet<float>("/cluster_histograms", HighFive::DataSpace::From(potential_energies_array));
    //     dataset.write(potential_energies_array);
    // }
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


struct CompareBySize{ template<typename T> bool operator()(const T& a, const T& b){return a.size()<b.size();}};

void DataCollector::try_cluster()
{
    clusters.setTarget(translator.particles.begin(),translator.particles.end());
    clusters.DBSCANrecursive(2, 1.4);


    {
        boost::multi_array<std::size_t,1> cluster_histogram(boost::extents[clusters.numClusters()]);
        for(std::size_t i = 0; i < clusters.numClusters(); ++i)
        {
            cluster_histogram[i] = (clusters.begin()+i)->size();
        }

        HighFive::DataSet dataset = FILE->createDataSet<std::size_t>("/cluster_histograms/time"+boost::lexical_cast<std::string>(timepoints.back()), HighFive::DataSpace::From(cluster_histogram));
        dataset.write(cluster_histogram);
    }

    tbb::concurrent_vector<ClusterVolumeParser> volumeParsers;
    tbb::parallel_for_each(clusters.begin(), clusters.end(), [&](const auto& cluster)
    {
        if(cluster.size() >= 20)
        {
            auto it = volumeParsers.emplace_back(cluster);
            it->setParameters(getParameters());
            it->parse();
        }
    });

    const float volume = PARALLEL_REDUCE(float, volumeParsers, [&](float i, const auto& parser)
    {
        return i + parser.result;
    });
    
    if(!volumeParsers.empty())
    std::max_element(std::begin(volumeParsers), std::end(volumeParsers), [](auto& a, auto& b)
    {
        return a.result < b.result;
    })->printXML("largest.vtu");
}