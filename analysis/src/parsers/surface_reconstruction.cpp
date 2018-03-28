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
    static const float extension = getParameters().analysis_cluster_volume_extension;
    check_for_aligned_box_setup();

    if(target_range->size() == 1)
    {
        volume = enhance::sphere_volume(1.f+extension);
        surface_area = enhance::sphere_surface(1.f+extension);
        return;
    }
    {
        // represent and manipulate 3D points
        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

        assert(target_range);
        // go over all subcluster members
        for( const auto& member : *target_range )
        {
            // generate a points
            cartesian point;
            point = member.position;

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
                // const cartesian point_to_check(Eigen::Vector3d(p.data()).cast<float>());
                // if(contains(point_to_check.template cast<float>()))
                    points->InsertNextPoint(p.data());
            }
        }
        
        // concrete dataset represents vertices, lines, polygons, and triangle strips
        vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
        polydata->SetPoints(points);

        vesDEBUG("vtkDelaunay3D")
        vtkSmartPointer<vtkDelaunay3D> delaunay = vtkSmartPointer<vtkDelaunay3D>::New();
        delaunay->SetInputData(polydata);
        delaunay->Update();

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