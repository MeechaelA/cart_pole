#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <time.h>
#include <vector>
#include "cart_pole.hpp"
#include "riccati_solver.hpp"
#include <cmath>
#include <iomanip>

class Simulation{
    public:
        Simulation();
        ~Simulation();
        bool start(std::string id, double end_time, unsigned int end_iteration, Eigen::MatrixXd Q, Eigen::MatrixXd R, Eigen::MatrixXd desired_state);

    private:
        unsigned int iteration;
        unsigned int end_iteration;


};

#endif
