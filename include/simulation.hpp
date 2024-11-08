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
        bool start(int trajectory_point, int simulation_num, double end_time, unsigned int end_iteration, Eigen::MatrixXd Q, Eigen::MatrixXd R, Eigen::MatrixXd start_state, Eigen::MatrixXd desired_state);
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
    std::tuple<std::vector<Eigen::MatrixXd>, std::vector<Eigen::MatrixXd>> create_power_2_trajectory(double initial_pos, int num_points, unsigned int dim_x, unsigned int dim_u);
    std::vector<double> double_range(double start, double end, int total);
    void output_simulation(std::string outfile, const std::vector<SimulationData>& m);
    void output_times(std::string outfile, int num_threads, int num_simulations, double total_time, std::vector<double> outer_times, std::vector<double> inner_times);
    void output_parent(std::string study_name);

}

#endif
