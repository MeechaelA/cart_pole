#include <iostream>
#include "simulation.hpp"
#include "cart_pole.hpp"

Simulation::Simulation(){
    this->iteration = 0;
    this->end_iteration = 0;
};

Simulation::~Simulation(){};

void Simulation::start(unsigned int end_iteration){
    this->end_iteration = end_iteration;

    while(iteration < end_iteration){
        std::cout << iteration << std::endl;
        iteration++;
    }
};