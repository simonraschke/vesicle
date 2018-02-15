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



ClusterVolumeParser::ClusterVolumeParser(const input_t& _cluster)
    : cluster(_cluster)
{
    result = 0;
    subclusters.setTarget(cluster.begin(), cluster.end());
    subclusters.DBSCANrecursive(1, 1.4);
}



void ClusterVolumeParser::parse()
{
    static const float extension = 1.4;

    for(const auto& subcluster : subclusters)
    {
        if(subcluster.size() < 6)
            continue;
        // represent and manipulate 3D points
        cartesian point;
        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

        vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
        sphereSource->SetRadius(1.f+extension);

        for( const auto& member_ptr : subcluster )
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
        vtkSmartPointer<vtkTriangleFilter> triangleFilter = vtkSmartPointer<vtkTriangleFilter>::New();
        triangleFilter->SetInputConnection(surfaceFilter->GetOutputPort());
            
        // estimate volume, area, shape index of triangle mesh
        vesLOG("vtkMassProperties")
        vtkSmartPointer<vtkMassProperties> massProperties = vtkSmartPointer<vtkMassProperties>::New();
        massProperties->SetInputConnection(triangleFilter->GetOutputPort());

        result += massProperties->GetVolume();
    }
}



// float ClusterVolumeParser::getVolume() const
// {
//     if(result)
//         return result->GetVolume();
//     else
//         throw std::logic_error("ClusterVolumeParser::structureProperties is nullptr");
// }



// float ClusterVolumeParser::getSurfaceArea() const
// {
//     if(result)
//         return result->GetSurfaceArea();
//     else
//         throw std::logic_error("ClusterVolumeParser::structureProperties is nullptr");
// }



// void ClusterVolumeParser::printXML(PATH system_complete_path) const
// {

//     if(!result)
//         throw std::logic_error("ClusterVolumeParser::triangleFilter is nullptr");

//     vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
//     writer->SetInputData(triangleFilter->GetOutput());
//     writer->SetFileName(system_complete_path.c_str());
//     writer->Write();
// }