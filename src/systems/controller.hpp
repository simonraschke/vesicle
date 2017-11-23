#pragma once

#include "systems/system.hpp"
#include "vesicleIO/history.hpp"
#include <tbb/parallel_invoke.h>
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
    static tbb::mutex signal_mutex;
};



struct SimulationControl
    : public Controller
{
    virtual void setup() override;
    virtual void start() override;
    virtual void pause() override;
    
private:
    HistoryStorage history_storage {};

    std::unique_ptr<tbb::flow::broadcast_node<tbb::flow::continue_msg>> start_node{nullptr};
    std::unique_ptr<tbb::flow::continue_node<tbb::flow::continue_msg>> step_node{nullptr};
    std::unique_ptr<tbb::flow::continue_node<tbb::flow::continue_msg>> thermostat_node{nullptr};
    std::unique_ptr<tbb::flow::continue_node<tbb::flow::continue_msg>> history_node{nullptr};
    std::unique_ptr<tbb::flow::continue_node<tbb::flow::continue_msg>> trajectory_node{nullptr};
};