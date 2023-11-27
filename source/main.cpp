#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <Eigen/Dense>
#include "simulation.hpp"
#include <vector>
#include <nlohmann/json.hpp>

// Command to start simulation
// ./cart_pole name 1 1 0.001 100.0 10
// [name num_threads num_simulations start_r end_r end_time]
int main(int argc, char *argv[]){

    // std::cout << "argc == " << argc << '\n';
    // for (int ndx{}; ndx != argc; ++ndx){
    //     std::cout << "argv[" << ndx << "] == " << std::quoted(argv[ndx]) << '\n';
    //     std::cout << "argv[" << argc << "] == " << static_cast<void*>(argv[argc]) << '\n';
    // }

    std::string study_name = argv[1];
    int num_threads = std::stoi(argv[2]);
    int num_simulations = std::stoi(argv[3]);
    double start_r = std::stod(argv[4]);
    double end_r = std::stod(argv[5]);
    double end_time = std::stod(argv[6]);
    int sims_ran = 0;

    omp_set_num_threads(num_threads);
    omp_set_nested(2);

    int end_iteration = 100000;
    std::vector<double> r_values = simulation_functions::double_range(start_r, end_r, num_simulations);

    unsigned int dim_x = 4;
    unsigned int dim_u = 1;
    std::vector<Eigen::MatrixXd> trajectory;
    std::vector<Eigen::MatrixXd> trajectory_prev;
    Eigen::MatrixXd start_state  = Eigen::MatrixXd::Zero(dim_x, dim_u);
    start_state(0) = 0.0;
    start_state(1) = 0.0;
    start_state(2) = 0.0;
    start_state(3) = 0.0;

    int num_pts = 30;
    double start = 0.0;
    double end = 500.0;

    // 
    // trajectory_prev = create_power_2_trajectory(start, num_pts, 4, 1);
    // trajectory = create_power_2_trajectory(trajectory_prev[1](0,0), num_pts, 4, 1);
    Eigen::MatrixXd point_zero  = Eigen::MatrixXd::Zero(dim_x, dim_u);
    point_zero(0) = 10.0;

    Eigen::MatrixXd point_one = Eigen::MatrixXd::Zero(dim_x, dim_u);
    point_one(0) = 30.0;

    Eigen::MatrixXd point_two  = Eigen::MatrixXd::Zero(dim_x, dim_u);
    point_two(0) = 50.0;

    Eigen::MatrixXd point_three  = Eigen::MatrixXd::Zero(dim_x, dim_u);
    point_three(0) = 70.0;

    std::tuple<std::vector<Eigen::MatrixXd>, std::vector<Eigen::MatrixXd>> value = simulation_functions::create_power_2_trajectory(0.0, 12, 4, 1);
    trajectory_prev = std::get<0>(value);
    trajectory = std::get<1>(value);

    omp_lock_t simulations_write_lock;
    omp_init_lock(&simulations_write_lock);
    std::vector<SimulationData> simulations;
 
    omp_lock_t local_simulations_write_lock;
    omp_init_lock(&local_simulations_write_lock);

    omp_lock_t statuses_write_lock;
    omp_init_lock(&statuses_write_lock);

    omp_lock_t local_statuses_write_lock;
    omp_init_lock(&local_statuses_write_lock);


    std::vector<double> inner_times;
    omp_lock_t inner_time_write_lock;
    omp_init_lock(&inner_time_write_lock);

    std::vector<double> outer_times;
    omp_lock_t outer_time_write_lock;
    omp_init_lock(&outer_time_write_lock);

    std::vector<TrajectoryStatus> statuses;
    std::cout << "Simulation Started" << std::endl;
    // Set Desired State & Pass to thread

    double total_time_initial = omp_get_wtime(); 
        #pragma omp parallel for num_threads(num_threads)
        for (int i = 0; i < trajectory.size(); i++)
        {
            std::vector<TrajectoryStatus> local_statuses;
            std::vector<SimulationData> local_simulations;

            double outer_time_initial = omp_get_wtime();

            Eigen::MatrixXd desired_state = trajectory[i];
            Eigen::MatrixXd prev_state = trajectory_prev[i];

            #pragma omp parallel for num_threads(num_threads) 
            for (int j = 0; j < num_simulations; j++)
            {
                double inner_time_initial = omp_get_wtime();
                
                bool success = false;
                Eigen::MatrixXd Q = Eigen::MatrixXd::Zero(dim_x, dim_x);
                Eigen::MatrixXd R = Eigen::MatrixXd::Zero(dim_u, dim_u);

                Q(0, 0) = 1.0;
                Q(1, 1) = 1.0;
                Q(2, 2) = 1.0;
                Q(3, 3) = 1.0;

                R(0, 0) = r_values[j];

                Simulation simulation = Simulation();
                success = simulation.start(i, j, end_time, end_iteration, Q, R, prev_state, desired_state);
                std::string status_message = std::to_string(i) + "," + std::to_string(j) + "," + std::to_string(success);

                TrajectoryStatus traj_status{
                    i,
                    j,
                    R(0,0),
                    simulation.status
                };

                omp_set_lock(&local_statuses_write_lock);
                local_statuses.push_back(traj_status);
                omp_unset_lock(&local_statuses_write_lock);

                omp_set_lock(&local_simulations_write_lock);
                local_simulations.push_back(simulation.get_data());
                omp_unset_lock(&local_simulations_write_lock);

                double inner_time_final = omp_get_wtime();
                double inner_time = inner_time_final - inner_time_initial;
                
                omp_set_lock(&inner_time_write_lock);
                inner_times.push_back(inner_time);
                omp_unset_lock(&inner_time_write_lock);
            }
            double outer_time_final = omp_get_wtime();
            double outer_time = outer_time_final - outer_time_initial;
            omp_set_lock(&inner_time_write_lock);
            outer_times.push_back(outer_time);
            omp_unset_lock(&inner_time_write_lock);

            for (int i_status = 0; i_status < local_statuses.size(); i_status++){
                omp_set_lock(&statuses_write_lock);
                statuses.push_back(local_statuses[i_status]);
                omp_unset_lock(&statuses_write_lock);                                    
            }

            for (int i_simulation = 0; i_simulation < local_simulations.size(); i_simulation++){
                omp_set_lock(&simulations_write_lock);
                simulations.push_back(local_simulations[i_simulation]);
                omp_unset_lock(&simulations_write_lock);                    
            }
        }
    
    double total_time_final = omp_get_wtime(); 
    double total_time_delta = total_time_final - total_time_initial;

    std::cout << "Total Simulation Time: " << total_time_delta << std::endl;
    std::cout << "Simulation Ended" << std::endl;
    std::cout << "Data Output Starting" << std::endl;

    simulation_functions::output_parent(study_name);
    simulation_functions::output_simulation(study_name, simulations);
    simulation_functions::output_times(study_name, num_threads, num_simulations, total_time_delta, outer_times, inner_times);

    std::cout << "Data Output Finished" << std::endl;

    return 0;
}