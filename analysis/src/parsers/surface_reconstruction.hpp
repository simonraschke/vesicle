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

#pragma once

#include "parser.hpp"
#include "cluster_parser.hpp"
#include <boost/filesystem.hpp>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkMassProperties.h>
#include <vtkSmartPointer.h>
#include <vtkDelaunay3D.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkTriangleFilter.h>
#include <vtkContourFilter.h>
#include <vtkSphereSource.h>
#include <vtkAppendFilter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>




class ClusterStructureParser
    : public Parser<vtkSmartPointer<vtkAppendFilter>>
{
public:
    typedef ClusterParser<PERIODIC::ON>::Cluster_t input_t;
    typedef input_t::MemberType::element_type member_t;
    typedef member_t::cartesian cartesian;
    typedef boost::filesystem::path PATH;

    explicit ClusterStructureParser(const input_t&);

    virtual void parse() override;

    float getVolume() const;
    float getSurfaceArea() const;

    void printXML(PATH) const;

protected:
    const input_t& cluster;
    float volume = 0;
    float surface_area = 0;
    ClusterParser<PERIODIC::OFF> subclusters;

private:
    // vtkSmartPointer<vtkAppendFilter> appendFilter {nullptr};
};