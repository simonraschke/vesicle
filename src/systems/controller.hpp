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



struct Controller
    : public ParameterDependentComponent
{
    virtual ~Controller() = default;

    virtual void setup() = 0;
    virtual void start() = 0;
    virtual void pause() = 0;
    
    static void signal(int SIG);

protected:
    Controller() = default;
    virtual void make_nodes() = 0;

    System system {};
    tbb::flow::graph flow {};

    static std::atomic<int> SIGNAL;
    static tbb::mutex signal_mutex;
};



struct SimulationControl
    : public Controller
{
    virtual void setup() override;
    virtual void start() override;
    virtual void pause() override;
    
protected:
    virtual void make_nodes() override;

private:
    HistoryStorage history_storage {};

    std::unique_ptr<tbb::flow::broadcast_node<tbb::flow::continue_msg>> start_node{nullptr};
    std::unique_ptr<tbb::flow::continue_node<tbb::flow::continue_msg>> step_node{nullptr};
    std::unique_ptr<tbb::flow::continue_node<tbb::flow::continue_msg>> thermostat_node{nullptr};
    std::unique_ptr<tbb::flow::continue_node<tbb::flow::continue_msg>> history_node{nullptr};
    std::unique_ptr<tbb::flow::continue_node<tbb::flow::continue_msg>> trajectory_node{nullptr};
};