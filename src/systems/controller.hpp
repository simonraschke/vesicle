#pragma once

#include "systems/system.hpp"
#include <tbb/flow_graph.h>
#include <csignal>



struct Controller
{
    virtual ~Controller() = default;

    virtual void setup() = 0;
    virtual void start() = 0;
    virtual void pause() = 0;
    
    static void signal(int SIG);

protected:
    Controller() = default;

    System system {};
    tbb::flow::graph flow {};

    static std::atomic<int> SIGNAL;
};



struct SimulationControl
    : public Controller
{
    virtual void setup() override;
    virtual void start() override;
    virtual void pause() override;

private:
    std::unique_ptr<tbb::flow::broadcast_node<tbb::flow::continue_msg>> start_node{nullptr};
    std::unique_ptr<tbb::flow::continue_node<tbb::flow::continue_msg>> step_node{nullptr};
    std::unique_ptr<tbb::flow::continue_node<tbb::flow::continue_msg>> trajectory_node{nullptr};
};