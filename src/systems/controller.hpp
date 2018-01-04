/*  
*   Copyright 2017 Simon Raschke
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

#include "systems/system.hpp"
#include "vesicleIO/history.hpp"
#include <tbb/parallel_invoke.h>
#include <tbb/flow_graph.h>
#include <csignal>
#include <chrono>



// a simulation control base class to derive from
// is parameter dependent component
// can catch signals if registered by std::signal 
struct Controller
    : public ParameterDependentComponent
{
    // destroy if derived is destroyed
    virtual ~Controller() = default;

    // derived MUST contain
    virtual void setup() = 0;
    virtual void start() = 0;
    virtual void pause() = 0;
    
    // static member function to catch signal
    // store in atomic which is accesible by derived
    static void signal(int SIG);

protected:
    // only to derive from
    Controller() = default;

    // force derived to somehow fill tbb::flow_graph
    virtual void make_nodes() = 0;

    // the actual system
    System system {};

    // the flow graph to fill with nodes 
    // stores the simulation pattern
    tbb::flow::graph flow {};

    // signal handling
    static std::atomic<int> SIGNAL;
    static tbb::mutex signal_mutex;
};




// derived from Controller
// this will control the simulation
// and implement virtual base class member functions
struct SimulationControl
    : public Controller
{
    // setup Controller::System and members
    virtual void setup() override;

    // start the flow graph and repeat until further notice
    virtual void start() override;

    // pause the simulation
    // TODO: implement
    virtual void pause() override;
    
protected:
    // fill the flow graph
    // by generating the node members
    virtual void make_nodes() override;

private:
    HistoryStorage history_storage {};

    std::unique_ptr<tbb::flow::broadcast_node<tbb::flow::continue_msg>> start_node{nullptr};
    std::unique_ptr<tbb::flow::continue_node<tbb::flow::continue_msg>> step_node{nullptr};
    std::unique_ptr<tbb::flow::continue_node<tbb::flow::continue_msg>> thermostat_node{nullptr};
    std::unique_ptr<tbb::flow::continue_node<tbb::flow::continue_msg>> history_node{nullptr};
    std::unique_ptr<tbb::flow::continue_node<tbb::flow::continue_msg>> trajectory_node{nullptr};
};