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
#include "enhance/concurrent_container.hpp"
#include "enhance/random.hpp"
#include "particles/particle_simple.hpp"
#include "geometries/grid.hpp"
#include <boost/filesystem.hpp>
#include <tbb/parallel_for_each.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
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
#include <vtkUnstructuredGrid.h>
#include <vtkVertexGlyphFilter.h>


// template<PERIODIC P>
// class ClusterBase;


// will construct a triangulated mesh from input
// if size of input points below threshold, a convex hull will be generated
//     this holds as a "good enough" estimation
// template<PERIODIC P>
class ClusterStructureParser
    : public Parser<vtkSmartPointer<vtkAppendFilter>>
{
public:
    typedef ParticleSimple member_t;
    typedef member_t::cartesian cartesian;
    typedef boost::filesystem::path PATH;

    // constructor
    // TODO: use setTarget
    ClusterStructureParser();
    
    void setTarget(enhance::ConcurrentDeque<ParticleSimple>&);

    // start the reconstruction
    virtual void parse() override;

    // output structure to PATH
    // call parse() beforehand
    void printXML(PATH) const;

    // get properties
    float getVolume() const;
    float getSurfaceArea() const;

protected:
    GridGeometry grid {};
    std::deque<bool> inside_flags {};
    float volume = 0;
    float surface_area = 0;
private:
    enhance::observer_ptr<enhance::ConcurrentDeque<ParticleSimple>> target_range {nullptr};

    // std::array<unsigned char,3> red = {255, 0, 0};
    // std::array<unsigned char,3> green = {0, 255, 0};
    // std::array<unsigned char,3> blue = {0, 0, 255};
    std::array<unsigned char,3> outside = {0, 0, 0};
    std::array<unsigned char,3> inside = {255, 153, 51};
};