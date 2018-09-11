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

#include "surface_reconstruction.hpp"


ClusterStructureParser::ClusterStructureParser()
{
    vtkObject::GlobalWarningDisplayOff();
    result = vtkSmartPointer<vtkAppendFilter>::New();
}



void ClusterStructureParser::setTarget(enhance::ConcurrentDeque<ParticleSimple>& target)
{
    vesDEBUG(__PRETTY_FUNCTION__)
    target_range = enhance::make_observer<enhance::ConcurrentDeque<ParticleSimple>>(&target);
}



void ClusterStructureParser::parse()
{
    check_for_aligned_box_setup();

    // tbb::task_arena limited(1);
    // limited.execute( [&]
    {
        const float points_per_sigma = 5;
        const float point_distance = getParameters().LJsigma / points_per_sigma;
        GridGeometry grid;
        {
            float x_center = 0.0;
            float y_center = 0.0;
            float z_center = 0.0;
            float x_edge = 0.0;
            float y_edge = 0.0;
            float z_edge = 0.0;
            {
                auto [min,max] = std::minmax_element(target_range->begin(), target_range->end(), [](const auto& p1, const auto& p2){ return p1.position(0) < p2.position(0); });
                x_center = (max->position(0) + min->position(0))/2;
                x_edge = max->position(0)-min->position(0) + getParameters().analysis_cluster_distance_threshold*2;
                // vesLOG(min->position(0) << " " << max->position(0));
                grid.x = x_edge * points_per_sigma + 1;
            }
            {
                auto [min,max] = std::minmax_element(target_range->begin(), target_range->end(), [](const auto& p1, const auto& p2){ return p1.position(1) < p2.position(1); });
                y_center = (max->position(1) + min->position(1))/2;
                y_edge = max->position(1)-min->position(1) + getParameters().analysis_cluster_distance_threshold*2;
                // vesLOG(min->position(1) << " " << max->position(1));
                grid.y = y_edge * points_per_sigma + 1;
            }
            {
                auto [min,max] = std::minmax_element(target_range->begin(), target_range->end(), [](const auto& p1, const auto& p2){ return p1.position(2) < p2.position(2); });
                z_center = (max->position(2) + min->position(2))/2;
                z_edge = max->position(2)-min->position(2) + getParameters().analysis_cluster_distance_threshold*2;
                // vesLOG(min->position(2) << " " << max->position(2));
                grid.z = z_edge * points_per_sigma + 1;
            }
            grid.generate();
            grid.scale(cartesian(point_distance, point_distance, point_distance));
            grid.shift(cartesian(x_center, y_center, z_center)-cartesian(x_edge, y_edge, z_edge)/2.f);
        }


        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
        for(const auto& member : *target_range)
        {
            points->InsertNextPoint(member.position.data());
        }
        auto polydata = vtkSmartPointer<vtkPolyData>::New();
        polydata->SetPoints(points);
        auto delaunay = vtkSmartPointer<vtkDelaunay3D>::New();
        delaunay->SetInputData(polydata);
        delaunay->Update();

        // append to overall mesh
        {
            // vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
            // polydata->SetPoints(points);
            result->AddInputData(polydata);
        }

        const auto squared_threshold = getParameters().analysis_cluster_distance_threshold*getParameters().analysis_cluster_distance_threshold;
        
        std::array<double, 3> point {};
        std::array<double, 3> pcoords {}; 
        std::array<double, 4> weights {};
        int subId;

        const auto inside_points = std::accumulate(std::begin(grid.points), std::end(grid.points), (std::size_t)0, [&](auto j, const auto& p)
        {   
            point = {p(0), p(1), p(2)};
            vtkIdType cellId = -1;
            if(target_range->size() > 6)
            {
                cellId = delaunay->GetOutput()->FindCell(point.data(), NULL, 0, .1, subId, pcoords.data(), weights.data());
            }

            if(cellId >= 0) 
                return j+1;
            else if(std::find_if(std::cbegin(*target_range), std::cend(*target_range), [&](const auto& mem){ return squared_distance(mem.position, p) <= squared_threshold; }) !=  std::end(*target_range))
                return j+1;
            else
                return j;
        });

        // const auto inside_points = tbb::parallel_reduce(tbb::blocked_range<typename decltype(grid.points)::const_iterator>( std::begin(grid.points), std::end(grid.points) ), (std::size_t)0 , [&](auto& r, std::size_t i) 
        // {

        //     return i + std::accumulate(std::begin(r), std::end(r), (std::size_t)0, [&](auto j, const auto& p)
        //     {   
        //         vtkIdType cellId = -1;
        //         if(target_range->size() > 6)
        //         {
        //             std::array<double, 3> point {p(0), p(1), p(2)};
        //             std::array<double, 3> pcoords {}; 
        //             std::array<double, 4> weights {};
        //             int subId;
        //             cellId = delaunay->GetOutput()->FindCell(point.data(), NULL, 0, .1, subId, pcoords.data(), weights.data());
        //         }
        //         // else

        //         if(cellId >= 0) 
        //             return j+1;
        //         else if(std::find_if(std::cbegin(*target_range), std::cend(*target_range), [&](const auto& mem){ return squared_distance(mem.position, p) <= squared_threshold; }) !=  std::end(*target_range))
        //             return j+1;
        //         else
        //             return j;
        //     });
        // }, std::plus<std::size_t>(), tbb::static_partitioner());

        volume += (float(inside_points) * std::pow(point_distance,3) );
    }
}



float ClusterStructureParser::getVolume() const
{
    if(result)
        return volume;
    else
        throw std::logic_error("ClusterStructureParser::result is nullptr. No structure");
}



float ClusterStructureParser::getSurfaceArea() const
{
    if(result)
        return surface_area;
    else
        throw std::logic_error("ClusterStructureParser::result is nullptr. No structure");
}



void ClusterStructureParser::printXML(PATH __attribute__((unused)) system_complete_path) const
{
    vesDEBUG(__PRETTY_FUNCTION__)

    if(!result)
        throw std::logic_error("ClusterStructureParser::result is nullptr. No structure");

    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetInputConnection(result->GetOutputPort());
    writer->SetFileName(system_complete_path.c_str());
    writer->Write();
}