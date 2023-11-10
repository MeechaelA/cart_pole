#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <Eigen/Dense>
#include "simulation.hpp"
#include <vector>

std::vector<double> double_range(int start, int end, double step){
    std::vector<double> values;
    for (int i=start; i < end; i++){
        values.push_back(i*step);
    }
    return values;
}

int main(int argc, char *argv[]){

    std::cout << "argc == " << argc << '\n';
 
    for (int ndx{}; ndx != argc; ++ndx){
        std::cout << "argv[" << ndx << "] == " << std::quoted(argv[ndx]) << '\n';
        std::cout << "argv[" << argc << "] == " << static_cast<void*>(argv[argc]) << '\n';
    }


    std::cout << "Simulation Started" << std::endl;
    int num_threads = std::stoi(argv[1]);
    omp_set_num_threads(num_threads);

    int end_iteration = 100000;
    double end_time = 10.0;

    std::vector<double> r_values = double_range(1, 100, 0.0001);
    int num_simulations = r_values.size();
    std::cout << "Num Sims: " << num_simulations << std::endl;

    #pragma omp parallel //num_threads(num_threads)
    {
        #pragma omp for //schedule(dynamic, 1)
            for (int i = 0; i < num_simulations; i++)
            {
                int tid = omp_get_thread_num();

                unsigned int dim_x = 4;
                unsigned int dim_u = 1;
                Eigen::MatrixXd Q = Eigen::MatrixXd::Zero(dim_x, dim_x);
                Eigen::MatrixXd R = Eigen::MatrixXd::Zero(dim_u, dim_u);
                Q(0, 0) = 1.0;
                Q(1, 1) = 1.0;
                Q(2, 2) = 1.0;
                Q(3, 3) = 10.0;
                R(0, 0) = r_values[i];

                Simulation simulation = Simulation();
                simulation.start(std::to_string(tid) + "_" + std::to_string(i), end_time, end_iteration, Q, R);
            }
    }

    std::cout << "Simulation Ended" << std::endl;



    return argc == 3 ? EXIT_SUCCESS : EXIT_FAILURE; // optional return value
}