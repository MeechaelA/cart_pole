#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <Eigen/Dense>
#include "simulation.hpp"
#include <vector>
#include <mpi.h>

std::vector<double> double_range(double start, double end, int total){
    std::vector<double> values;
    double step = (end-start)/total;
    for (int i=1; i <= total; i++){
        values.push_back(i*step);
    }
    return values;
}

struct TrajectoryStatus{
    int processor_id = -1;
    int thread_id = -1;
    int trajectory_id = -1;
    bool trajectory_status = 0; // Some day could be a range of values to report "Goodness..." (Do when transition to rust...)
};

int main(int argc, char *argv[]){

    // mpirun -np 4 ./cart_pole 10 10 0.001 100 25.0


    int numprocs,rank;

    // Command to start simulation
    // 
    // ./cart_pole 10 10 0.001 100.0
    // [num_threads num_simulations start_r end_r end_time]

    int num_threads = std::stoi(argv[1]);
    int num_simulations = std::stoi(argv[2]);
    double start_r = std::stod(argv[3]);
    double end_r = std::stod(argv[4]);
    double end_time = std::stod(argv[5]);

    omp_set_num_threads(num_threads);

    int end_iteration = 100000;
    std::vector<double> r_values = double_range(start_r, end_r, num_simulations);

    unsigned int dim_x = 4;
    unsigned int dim_u = 1;

    std::vector<Eigen::MatrixXd> trajectory;
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
       
    std::vector<TrajectoryStatus> statuses; 


    MPI_Init(&argc,&argv); // Initalize MPI environment
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs); //get total number of processes
    MPI_Comm_rank(MPI_COMM_WORLD,&rank); // get process identity number
    // MPI_Buffer_attach(statuses, trajectory.size() * num_simulations);

    int num_per_rank = std::floor((trajectory.size()/numprocs));
    if (trajectory.size()%numprocs > 0){
        num_per_rank += 1;
    }


    MPI_Status status;


    // Make special type for custom message passing...
    MPI_Datatype mpi_trajectory_message;
    TrajectoryStatus trajectory_status;
    int B[] = {1,1,1,1};
    MPI::Aint D[] = {
        offsetof(struct TrajectoryStatus, processor_id),
        offsetof(struct TrajectoryStatus, thread_id),
        offsetof(struct TrajectoryStatus, trajectory_id),
        offsetof(struct TrajectoryStatus, trajectory_status),
        sizeof(struct TrajectoryStatus)
    };
    MPI_Datatype T[4] = {
        MPI::INT,
        MPI::INT,
        MPI::INT,
        MPI::BOOL
    };
    MPI_Type_create_struct(4, B, D, T, &mpi_trajectory_message);
    MPI_Type_commit(&mpi_trajectory_message);

    // 
    MPI_Request request = MPI_REQUEST_NULL;

    std::cout << "argc == " << argc << '\n';
    for (int ndx{}; ndx != argc; ++ndx){
        std::cout << "argv[" << ndx << "] == " << std::quoted(argv[ndx]) << '\n';
        std::cout << "argv[" << argc << "] == " << static_cast<void*>(argv[argc]) << '\n';
    }


    std::cout << "Simulation Started" << std::endl;

    // Set Desired State & Pass to thread
    for (int i = rank * num_per_rank; i <= rank*num_per_rank + num_per_rank; i++)
    {
        TrajectoryStatus local_statuses[num_simulations];
        if (i < trajectory.size()){
            std::cout << i << "," << num_per_rank << std::endl;

            #pragma omp parallel
            {
                int tid =  omp_get_thread_num();
                // Need to tell this thread if it was successful or not for given time to
                Eigen::MatrixXd desired_state = trajectory[i];
                
                #pragma omp for //schedule(dynamic, 1)
                for (int j = 0; j < num_simulations; j++){

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
                    success = simulation.start(std::to_string(i) + "_" + std::to_string(j), end_time, end_iteration, Q, R, desired_state);
                    std::string status_message = std::to_string(i) + "," + std::to_string(j) + "," + std::to_string(success);

                    TrajectoryStatus traj_status{
                        rank,
                        tid,
                        j,
                        success
                    };
                    local_statuses[j] = traj_status;
                }
                // #TrajectoryStatus recv_status;
                // #MPI_Recv(&recv_status, 1, mpi_trajectory_message, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                // #std::cout << "Receive: " << std::to_string(recv_status.processor_id) << "," << std::to_string(recv_status.thread_id) << "," << std::to_string(recv_status.trajectory_id) << "," << std::boolalpha << recv_status.trajectory_status << std::endl;
                // #statuses[recv_status.processor_id * recv_status.trajectory_id] = recv_status;
            }
            MPI_Isend(local_statuses, num_simulations, mpi_trajectory_message, 0, 0, MPI_COMM_WORLD, &request); 
            MPI_Irecv(local_statuses, num_simulations, mpi_trajectory_message, rank, 1, MPI_COMM_WORLD, &request);
            for (int i_status = 0; i < sizeof(local_statuses)/sizeof(local_statuses[0]); i++){
                statuses.push_back(local_statuses[i]);
            }
        }
        else{
            std::cout << i << "," << "NO OP" << std::endl;
        }


    }
    std::cout << "Simulation Ended" << std::endl;




    // Trajectory Point #,  Simulation #, Status
    std::ofstream outfile;
    outfile.open(std::to_string(rank) + "simulation.csv");
    for (int i = 0; i < statuses.size(); i++){
        outfile << std::to_string(statuses[i].processor_id) << "," << std::to_string(statuses[i].thread_id) << "," << std::to_string(statuses[i].trajectory_id) << "," << std::boolalpha << statuses[i].trajectory_status << std::endl;
    }

    MPI_Finalize();   

    return 0;
}