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
#include "vtkPowerCrustSurfaceReconstruction.h"




// will construct a triangulated mesh from input
// if size of input points below threshold, a convex hull will be generated
//     this holds as a "good enough" estimation
// WARNING: if size beyond threshold, powercrust will generate the mesh
//     serious memory leak
//     extremely good results
//     submit job with loads of memory
class ClusterStructureParser
    : public Parser<vtkSmartPointer<vtkAppendFilter>>
{
public:
    typedef ClusterParser<PERIODIC::ON>::Cluster_t input_t;
    typedef input_t::member_t::element_type member_t;
    typedef member_t::cartesian cartesian;
    typedef boost::filesystem::path PATH;

    // constructor
    // TODO: use setTarget
    explicit ClusterStructureParser(const input_t&);

    // start the reconstruction
    virtual void parse() override;

    // output structure to PATH
    // call parse() beforehand
    void printXML(PATH) const;

    // get properties
    float getVolume() const;
    float getSurfaceArea() const;
    float getOrder() const;
    cartesian getCenter() const;
    std::size_t getNumMembers() const;

    // member evaluation
    template<PARTICLETYPE P>
    bool containsMemberType() const;
    template<PARTICLETYPE P>
    bool notContainsMemberType() const;

protected:
    const input_t& cluster;
    float volume = 0;
    float surface_area = 0;
    ClusterParser<PERIODIC::OFF> subclusters;

private:
    std::size_t powercrust_min = 110; // TODO: make it input parameter
    // vtkSmartPointer<vtkAppendFilter> appendFilter {nullptr};
};



template<PARTICLETYPE P>
bool ClusterStructureParser::containsMemberType() const
{
    return std::any_of(std::begin(cluster), std::end(cluster),[](const auto& p){ return P == p->type; });
}



template<PARTICLETYPE P>
bool ClusterStructureParser::notContainsMemberType() const
{
    return std::none_of(std::begin(cluster), std::end(cluster),[](const auto& p){ return P == p->type; });
}