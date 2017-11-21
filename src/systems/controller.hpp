#pragma once

#include <tbb/flow_graph.h>



struct Controller
{
    virtual ~Controller() = default;

    virtual void start() = 0;
    virtual void pause() = 0;

protected:
    Controller() = default;

    tbb::flow::graph flow;
};



struct SimulationControl
    : public Controller
{
    virtual void start() override;
    virtual void pause() override;
};