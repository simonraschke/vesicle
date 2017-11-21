#include <cstdlib>
#include <iostream>
#include "systems/controller.hpp"


int main()
{
    SimulationControl control;
    control.setup();
    control.start();

    return EXIT_SUCCESS;
}
