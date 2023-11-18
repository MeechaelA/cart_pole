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
#include <omp.h>
#include <map>

struct TrajectoryStatus{
    int trajectory_id;
    int simulation_id;
    double r_value;
    double trajectory_status; // Some day could be a range of values to report "Goodness..." (Do when transition to rust...)
};

struct SimulationData{
    TrajectoryStatus status;
    std::map<std::string, std::map<std::string, std::vector<double>>> data;
};

class Simulation{
    public:
        Simulation();
        ~Simulation();
        bool start(std::string id, double end_time, unsigned int end_iteration, Eigen::MatrixXd Q, Eigen::MatrixXd R, Eigen::MatrixXd start_state, Eigen::MatrixXd desired_state);
        SimulationData get_data();
        std::string id;
        double status;


        
    private:
        unsigned int iteration;
        unsigned int end_iteration;

        SimulationData simulation_data;
        std::vector<double> times;
        std::vector<double> cart_positions;
        std::vector<double> cart_velocities;
        std::vector<double> cart_accelerations;
        std::vector<double> pole_positions;
        std::vector<double> pole_velocities;
        std::vector<double> pole_accelerations;
        std::vector<double> forces;
};

namespace simulation_functions{
    void output(std::string outfile, const std::vector<SimulationData>& m);
}

#endif
