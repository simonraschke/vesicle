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



ClusterStructureParser::ClusterStructureParser(const input_t& _cluster)
    : cluster(_cluster)
    // , appendFilter(vtkSmartPointer<vtkAppendFilter>::New())
{
    vtkObject::GlobalWarningDisplayOff();
    result = vtkSmartPointer<vtkAppendFilter>::New();
    // result = 0;
    subclusters.setTarget(cluster.begin(), cluster.end());
}



void ClusterStructureParser::parse()
{
    if(getParameters().analysis_cluster_algorithm == "DBSCAN")
        subclusters.DBSCANrecursive(getParameters().analysis_cluster_minimum_size, getParameters().analysis_cluster_distance_threshold);
    else
        vesCRITICAL("cluster algorithm " << getParameters().analysis_cluster_algorithm << " unknown")
    
    static const float extension = getParameters().analysis_cluster_volume_extension;
    check_for_aligned_box_setup();

    // go over all subclusters
    for(const auto& subcluster : subclusters)
    {
        if(subcluster.size() < 3)
            continue;

        // represent and manipulate 3D points
        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

        // go over all subcluster members
        for( const auto& member_ptr : subcluster )
        {
            // generate a points
            cartesian point;
            point = member_ptr->position;

            // generate a sphere and place it around point
            vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
            sphereSource->SetRadius(1.f+extension);
            sphereSource->SetCenter( point(0), point(1), point(2) );
            sphereSource->SetThetaResolution(10);
            sphereSource->SetPhiResolution(6);
            sphereSource->Update();

            // save the new positioned sphere in vtkPolyData
            vtkSmartPointer<vtkPolyData> spherePolyData = sphereSource->GetOutput();

            // store all points in vtkPoints
            for(vtkIdType i = 0; i < spherePolyData->GetNumberOfPoints(); i++)
            {
                std::array<double,3> p{};

                // store point values in array
                spherePolyData->GetPoint(i,p.data());

                // insert into points container if inside of box
                // const Eigen::Vector3d point_to_check(p.data());
                const cartesian point_to_check(Eigen::Vector3d(p.data()).cast<float>());
                if(contains(point_to_check.cast<float>()))
                    points->InsertNextPoint(p.data());
            }
        }
        
        // concrete dataset represents vertices, lines, polygons, and triangle strips
        vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
        polydata->SetPoints(points);

        // convex hull
        vesDEBUG("vtkDelaunay3D")
        vtkSmartPointer<vtkDelaunay3D> delaunay = vtkSmartPointer<vtkDelaunay3D>::New();
        delaunay->SetInputData(polydata);

        // extract outer (polygonal) surface
        vesDEBUG("vtkDataSetSurfaceFilter")
        vtkSmartPointer<vtkDataSetSurfaceFilter> surfaceFilter = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
        surfaceFilter->SetInputConnection(delaunay->GetOutputPort());

        // extract outer (polygonal) surface
        vesDEBUG("vtkTriangleFilter")
        vtkSmartPointer<vtkTriangleFilter> triangleFilter = vtkSmartPointer<vtkTriangleFilter>::New();
        triangleFilter->SetInputConnection(surfaceFilter->GetOutputPort());
            
        // append to overall mesh
        result->AddInputConnection(triangleFilter->GetOutputPort());
        
        // estimate volume, area, shape index of triangle mesh
        vesDEBUG("vtkMassProperties")
        vtkSmartPointer<vtkMassProperties> massProperties = vtkSmartPointer<vtkMassProperties>::New();
        massProperties->SetInputConnection(surfaceFilter->GetOutputPort());
        massProperties->Update();

        volume += massProperties->GetVolume();
        surface_area += massProperties->GetSurfaceArea();
    }
    result->Update();
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



void ClusterStructureParser::printXML(PATH system_complete_path) const
{
    vesDEBUG(__PRETTY_FUNCTION__)

    if(!result)
        throw std::logic_error("ClusterStructureParser::result is nullptr. No structure");

    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetInputConnection(result->GetOutputPort());
    writer->SetFileName(system_complete_path.c_str());
    writer->Write();
}



ClusterStructureParser::cartesian ClusterStructureParser::getCenter() const
{
    return std::accumulate(std::begin(cluster), std::end(cluster), cartesian(cartesian::Zero()), [](const cartesian& c, const auto& p){ return c + p->position; }) / cluster.size();
}



std::size_t ClusterStructureParser::getNumMembers() const
{
    return cluster.size();
}



float ClusterStructureParser::getOrder() const
{
    const auto center_ = getCenter();
    return std::accumulate(std::begin(cluster), std::end(cluster), float(0), [&](float order, const auto& p)
    { 
        return order + distanceVector(center_,p->position).normalized().dot(p->orientation.normalized());
    }) / cluster.size();
}