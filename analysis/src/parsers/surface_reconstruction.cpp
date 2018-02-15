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



SurfaceReconstructionParser::SurfaceReconstructionParser(const input_t& _cluster)
    : cluster(_cluster)
{
    
}



void SurfaceReconstructionParser::parse()
{
    
    // represent and manipulate 3D points
    cartesian point;
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

    static const float extension = 1.4;
    vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
    sphereSource->SetRadius(1.f+extension);

    for( const auto& member_ptr : cluster )
    {
        point = member_ptr->position;
        sphereSource->SetCenter( point(0), point(1), point(2) );
        sphereSource->Update();

        vtkSmartPointer<vtkPolyData> spherePolyData = sphereSource->GetOutput();

        std::array<double,3> p{};
        for(vtkIdType i = 0; i < spherePolyData->GetNumberOfPoints(); i++)
        {
            spherePolyData->GetPoint(i,p.data());
            points->InsertNextPoint(p.data());
        }
    }
    
    // concrete dataset represents vertices, lines, polygons, and triangle strips
    vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
    polydata->SetPoints(points);

    // convex hull
    vesLOG("vtkDelaunay3D")
    vtkSmartPointer<vtkDelaunay3D> delaunay = vtkSmartPointer<vtkDelaunay3D>::New();
    delaunay->SetInputData(polydata);

    // extract outer (polygonal) surface
    vesLOG("vtkDataSetSurfaceFilter")
    vtkSmartPointer<vtkDataSetSurfaceFilter> surfaceFilter = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
    surfaceFilter->SetInputConnection(delaunay->GetOutputPort());

    // extract outer (polygonal) surface
    vesLOG("vtkTriangleFilter")
    triangleFilter = vtkSmartPointer<vtkTriangleFilter>::New();
    triangleFilter->SetInputConnection(surfaceFilter->GetOutputPort());
        
    // estimate volume, area, shape index of triangle mesh
    vesLOG("vtkMassProperties")
    result = vtkSmartPointer<vtkMassProperties>::New();
    result->SetInputConnection(triangleFilter->GetOutputPort());
}



float SurfaceReconstructionParser::getVolume() const
{
    if(result)
        return result->GetVolume();
    else
        throw std::logic_error("SurfaceReconstructionParser::structureProperties is nullptr");
}



float SurfaceReconstructionParser::getSurfaceArea() const
{
    if(result)
        return result->GetSurfaceArea();
    else
        throw std::logic_error("SurfaceReconstructionParser::structureProperties is nullptr");
}



void SurfaceReconstructionParser::printXML(PATH system_complete_path) const
{

    if(!result)
        throw std::logic_error("SurfaceReconstructionParser::triangleFilter is nullptr");

    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetInputData(triangleFilter->GetOutput());
    writer->SetFileName(system_complete_path.c_str());
    writer->Write();
}