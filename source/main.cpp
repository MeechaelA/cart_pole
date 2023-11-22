#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <Eigen/Dense>
#include "simulation.hpp"
#include <vector>
#include <nlohmann/json.hpp>


std::vector<Eigen::MatrixXd> create_trajectory(double initial_pos, double end_position, int num_points, unsigned int dim_x, unsigned int dim_u){
    std::vector<Eigen::MatrixXd> trajectory;
    Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(dim_x, dim_u);
    for (int i = 0; i < num_points; i++){
        matrix(0) += (end_position - initial_pos) / num_points;
        trajectory.push_back(matrix);
    }
    return trajectory;
}


std::vector<double> double_range(double start, double end, int total){
    std::vector<double> values;
    double step = (end-start)/total;
    for (int i=1; i <= total; i++){
        values.push_back(i*step);
    }
    return values;
}

// mpirun -np 4 ./cart_pole 10 10 0.001 100 25.0
// Command to start simulation
//
// ./cart_pole 10 10 0.001 100.0
// [num_threads num_simulations start_r end_r end_time]
int main(int argc, char *argv[]){

    // std::cout << "argc == " << argc << '\n';
    // for (int ndx{}; ndx != argc; ++ndx){
    //     std::cout << "argv[" << ndx << "] == " << std::quoted(argv[ndx]) << '\n';
    //     std::cout << "argv[" << argc << "] == " << static_cast<void*>(argv[argc]) << '\n';
    // }

    int num_threads = std::stoi(argv[1]);
    int num_simulations = std::stoi(argv[2]);
    double start_r = std::stod(argv[3]);
    double end_r = std::stod(argv[4]);
    double end_time = std::stod(argv[5]);
    int sims_ran = 0;

    omp_set_num_threads(num_threads);
    omp_set_nested(1);

    int end_iteration = 100000;
    std::vector<double> r_values = double_range(start_r, end_r, num_simulations);

    unsigned int dim_x = 4;
    unsigned int dim_u = 1;
    std::vector<Eigen::MatrixXd> trajectory;
    std::vector<Eigen::MatrixXd> trajectory_prev;
    Eigen::MatrixXd start_state  = Eigen::MatrixXd::Zero(dim_x, dim_u);

    start_state(0) = 0.0;
    start_state(1) = 0.0;
    start_state(2) = 0.0;
    start_state(3) = 0.0;

    int num_pts = 10;
    double start = 0.0;
    double end = 100.0;
    trajectory = create_trajectory(start + end/num_pts, end, num_pts, 4, 1);
    trajectory_prev = create_trajectory(start, end, num_pts, 4, 1);

    // trajectory.push_back(point_four);
    // trajectory.push_back(point_five);
    // trajectory.push_back(point_six);
    // trajectory.push_back(point_seven);
    // trajectory.push_back(point_eight);
    // trajectory.push_back(point_nine);

    omp_lock_t simulations_write_lock;
    omp_init_lock(&simulations_write_lock);
    std::vector<SimulationData> simulations;
 
    omp_lock_t local_simulations_write_lock;
    omp_init_lock(&local_simulations_write_lock);

    omp_lock_t statuses_write_lock;
    omp_init_lock(&statuses_write_lock);

    omp_lock_t local_statuses_write_lock;
    omp_init_lock(&local_statuses_write_lock);

    std::vector<TrajectoryStatus> statuses;
    std::cout << "Simulation Started" << std::endl;
    // Set Desired State & Pass to thread
    #pragma omp parallel
    {
        std::vector<TrajectoryStatus> local_statuses;
        std::vector<SimulationData> local_simulations;
        #pragma omp for
        for (int i = 0; i < trajectory.size(); i++)
        {

            Eigen::MatrixXd desired_state = trajectory[i];
            Eigen::MatrixXd prev_state = trajectory_prev[i];

            #pragma omp parallel
            {
                #pragma omp for
                for (int j = 0; j < num_simulations; j++)
                {

                    
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
                }
            }
        }
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

    std::cout << "Simulation Ended" << std::endl;

    std::cout << "Data Output Starting" << std::endl;
    simulation_functions::output("simulation.json", simulations);  
    std::cout << "Data Output Finished" << std::endl;

    return 0;
}