#include "algorithm.hpp"



Algorithm::Algorithm ()
    : target_range({}) 
{

}



void Algorithm::setTarget(std::initializer_list<target_type> list)
{
    target_range = list;
}



// void Algorithm::setParameters(Parameters prms)
// {
//     parameters = std::make_unique<Parameters>(prms);
// }