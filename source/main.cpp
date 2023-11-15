#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <Eigen/Dense>
#include "simulation.hpp"
#include <vector>

std::vector<double> double_range(double start, double end, int total){
    std::vector<double> values;
    double step = (end-start)/total;
    for (int i=1; i <= total; i++){
        values.push_back(i*step);
    }
    return values;
}

struct TrajectoryStatus{
    int trajectory_id = -1;
    int simulation_id = -1;
    double r_value = 0.0;
    bool trajectory_status = 0; // Some day could be a range of values to report "Goodness..." (Do when transition to rust...)
};

int main(int argc, char *argv[]){

    // mpirun -np 4 ./cart_pole 10 10 0.001 100 25.0
    // Command to start simulation
    // 
    // ./cart_pole 10 10 0.001 100.0
    // [num_threads num_simulations start_r end_r end_time]




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
    Eigen::MatrixXd start_state  = Eigen::MatrixXd::Zero(dim_x, dim_u);
    start_state(0) = 0.0;
    start_state(1) = 0.0;
    start_state(2) = 3.14159;
    start_state(3) = 0.0;
    Eigen::MatrixXd prev_state = start_state;

    Eigen::MatrixXd point_zero  = Eigen::MatrixXd::Zero(dim_x, dim_u);
    point_zero(0) = 0.0;

    Eigen::MatrixXd point_one   = Eigen::MatrixXd::Zero(dim_x, dim_u);
    point_one(0) = 500.0;

    Eigen::MatrixXd point_two   = Eigen::MatrixXd::Zero(dim_x, dim_u);
    point_two(0) = -500.0;

    Eigen::MatrixXd point_three = Eigen::MatrixXd::Zero(dim_x, dim_u);
    point_three(0) = 0.0;

    Eigen::MatrixXd point_four  = Eigen::MatrixXd::Zero(dim_x, dim_u);
    point_four(0) = 1000.0;

    Eigen::MatrixXd point_five  = Eigen::MatrixXd::Zero(dim_x, dim_u);
    point_five(0) = -1000.0;

    Eigen::MatrixXd point_six   = Eigen::MatrixXd::Zero(dim_x, dim_u);
    point_six(0) = 0.0;
    
    trajectory.push_back(point_zero);
    trajectory.push_back(point_one); 
    trajectory.push_back(point_two);
    trajectory.push_back(point_three);
    trajectory.push_back(point_four);
    trajectory.push_back(point_five);
    trajectory.push_back(point_six);
       
    omp_lock_t write_lock;
    omp_init_lock(&write_lock);
    std::vector<TrajectoryStatus> statuses; 

    std::cout << "argc == " << argc << '\n';
    for (int ndx{}; ndx != argc; ++ndx){
        std::cout << "argv[" << ndx << "] == " << std::quoted(argv[ndx]) << '\n';
        std::cout << "argv[" << argc << "] == " << static_cast<void*>(argv[argc]) << '\n';
    }


    std::cout << "Simulation Started" << std::endl;
    // Set Desired State & Pass to thread
    #pragma omp parallel
    {
        #pragma omp for
        for (int i = 0; i <= trajectory.size(); i++)
        {
            TrajectoryStatus local_statuses[num_simulations];
            if (i < trajectory.size()){
                sims_ran++;
                // std::cout << i << "," << num_per_rank << std::endl;
                // std::vector<TrajectoryStatus> local_statuses;
                // Need to tell this thread if it was successful or not for given time to
                if (i > 0){
                    prev_state = trajectory[i-1];
                }
                Eigen::MatrixXd desired_state = trajectory[i];
                int tid =  omp_get_thread_num();
                
                
                #pragma omp parallel
                {
                    #pragma omp for //schedule(dynamic, 1)
                    for (int j = 0; j < num_simulations; j++)
                    {
                        bool success = false;
                        int tid = omp_get_thread_num();
                        Eigen::MatrixXd Q = Eigen::MatrixXd::Zero(dim_x, dim_x);
                        Eigen::MatrixXd R = Eigen::MatrixXd::Zero(dim_u, dim_u);
                        Q(0, 0) = 1.0;
                        Q(1, 1) = 1.0;
                        Q(2, 2) = 1.0;
                        Q(3, 3) = 1.0;
                        R(0, 0) = r_values[j];

                        Simulation simulation = Simulation();
                        success = simulation.start(std::to_string(i) + "_" + std::to_string(j), end_time, end_iteration, Q, R, prev_state, desired_state);
                        std::string status_message = std::to_string(i) + "," + std::to_string(j) + "," + std::to_string(success);

                        TrajectoryStatus traj_status{
                            i,
                            j,
                            R(0,0),
                            success
                        };
                        omp_set_lock(&write_lock);
                        statuses.push_back(traj_status);
                        omp_unset_lock(&write_lock);
                    }
                }
            }
            else{
                std::cout << i << "," << "NO OP" << std::endl;
            }
        }
    }
    omp_destroy_lock(&write_lock);
    std::cout << "Simulation Ended" << std::endl;

    // Trajectory Point #,  Simulation #, Status
    std::ofstream outfile;
    outfile.open("simulation.csv");
    outfile << "trajectory" << "," << "simulation" << "," << "input_cost" << "," << "status" << std::endl;
    for (int i = 0; i < statuses.size(); i++){
        outfile << std::to_string(statuses[i].trajectory_id) << "," << std::to_string(statuses[i].simulation_id) << "," << std::to_string(statuses[i].r_value) << "," << std::boolalpha << statuses[i].trajectory_status << std::endl;
    }

    // MPI_Finalize();   

    return 0;
}