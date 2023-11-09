
#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <time.h>
#include <vector>
#include "simulation.hpp"
#include "cart_pole.hpp"
#include "riccati_solver.hpp"
#include <cmath>
#include <iomanip>

#define PRINT_MAT(X) std::cout << #X << ":\n" << X << std::endl << std::endl

Simulation::Simulation(){
    this->iteration = 0;
    this->end_iteration = 0;
};

Simulation::~Simulation(){};

void Simulation::start(unsigned int end_iteration){
    std::ofstream outfile;
    outfile.open("output.csv");

    this->end_iteration = end_iteration;



    const uint dim_x = 4;
    const uint dim_u = 1;

    Eigen::MatrixXd X = Eigen::MatrixXd::Zero(dim_x, dim_u);
    Eigen::MatrixXd X_prev = Eigen::MatrixXd::Zero(dim_x, dim_u);
    Eigen::MatrixXd X_dot = Eigen::MatrixXd::Zero(dim_x, dim_u);
    Eigen::MatrixXd X_dot_prev = Eigen::MatrixXd::Zero(dim_x, dim_u);
    Eigen::MatrixXd X_desired = Eigen::MatrixXd::Zero(dim_x, dim_u);
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(dim_x, dim_x);
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(dim_x, dim_u);
    Eigen::MatrixXd Q = Eigen::MatrixXd::Zero(dim_x, dim_x);
    Eigen::MatrixXd R = Eigen::MatrixXd::Zero(dim_u, dim_u);
    Eigen::MatrixXd P = Eigen::MatrixXd::Zero(dim_x, dim_x);
    Eigen::MatrixXd K_continous = Eigen::MatrixXd::Zero(dim_x, dim_x);
    Eigen::MatrixXd K_discrete = Eigen::MatrixXd::Zero(dim_x, dim_x);
    Eigen::MatrixXd u = Eigen::MatrixXd::Zero(dim_u, dim_u);

    X(0) = 0.0;
    X(1) = 0.0;
    X(2) = 0.0;
    X(3) = 0.0;
    X_prev = X;

    X_dot(0) = 0.0;
    X_dot(1) = 0.0;
    X_dot(2) = 0.0;
    X_dot(3) = 0.0;
    X_dot_prev = X_dot;

    X_desired(0) = 0.0;
    X_desired(1) = 0.0;
    X_desired(2) = 3.14159/2.0;
    X_desired(3) = 0.0;



    CartPole    cart_pole;

    cart_pole.friction = 1.0;
    cart_pole.cart_mass = 5.0;
    cart_pole.pole_mass = 1.0;
    cart_pole.pole_length = 2.0;
    double b = 1.0;
    

    double pole_dis = 3.14159/2.0 + 0.01;
    double cart_dis = 0.0;
    double pole_vel = 0.0;
    double cart_vel = 0.0;
    double pole_accel = 0.0;
    double cart_accel = 0.0;

    double pole_dis_prev = pole_dis;
    double cart_dis_prev = cart_dis;
    double pole_vel_prev = pole_vel;
    double cart_vel_prev = cart_vel;
    double pole_accel_prev = pole_accel;
    double cart_accel_prev = cart_accel;

    double force = 0.0;
    double time = 0.0;
    double time_delta = 0.0001;
    double time_prev = 0.0;

    outfile.precision(12);
    outfile <<  "time" << "," << "cart_accel" << "," << "cart_vel" << "," << "cart_dis" << "," << "pole_accel" << "," << "pole_vel" << "," << "pole_dis" << "," <<  "force" << std::endl;

    while(iteration < end_iteration){
        A(0, 0) = 0.0;
        A(0, 1) = 1.0;
        A(0, 2) = 0.0;
        A(0, 3) = 0.0;

        A(1, 0) = 0.0;
        A(1, 1) = -1.0 * cart_pole.friction / cart_pole.cart_mass;
        A(1, 2) = b * cart_pole.pole_mass * cart_pole.gravity / cart_pole.cart_mass;
        A(1, 3) = 0.0;

        A(2, 0) = 0.0;
        A(2, 1) = 0.0;
        A(2, 2) = 0.0;
        A(2, 3) = 1.0;

        A(3, 0) = 0.0;
        A(3, 1) = -b * cart_pole.friction / cart_pole.cart_mass * cart_pole.pole_length;
        A(3, 2) = -b * (cart_pole.pole_mass + cart_pole.cart_mass)*cart_pole.gravity / cart_pole.cart_mass * cart_pole.pole_length;
        A(3, 3) = 1.0;

        B(0, 0) = 1.0;
        B(1, 0) = 1.0 / cart_pole.cart_mass;
        B(2, 0) = b / cart_pole.cart_mass*cart_pole.pole_length;
        B(3, 0) = 1.0;

        Q(0, 0) = 1.0;
        Q(1, 1) = 1.0;
        Q(2, 2) = 10.0;
        Q(3, 3) = 100.0;

        R(0, 0) = 0.001;

        double pole_dis_cos = cos(pole_dis);
        double pole_dis_sin = sin(pole_dis);
        // std::cout << "cos: "<< pole_dis_cos << std::endl;
        // std::cout << "sin: "<< pole_dis_cos << std::endl;
        double p = cart_pole.pole_mass * cart_pole.pole_length * cart_pole.pole_length * (cart_pole.cart_mass + cart_pole.pole_mass*pole_dis_cos);

        cart_accel = (1.0/p)*(pow(-cart_pole.pole_mass,2.0)*pow(cart_pole.pole_length,2.0)*cart_pole.gravity*pole_dis_cos*pole_dis_sin + cart_pole.pole_mass*pow(cart_pole.pole_length,2.0)*(cart_pole.pole_mass*cart_pole.pole_length*pow(pole_vel,2.0)*pole_dis_sin - cart_pole.friction*cart_vel)) + cart_pole.pole_mass*cart_pole.pole_length*cart_pole.pole_length*(1.0/p)*force;
        pole_accel = (1.0/p)*((cart_pole.pole_mass+cart_pole.cart_mass)*cart_pole.pole_mass*cart_pole.gravity*cart_pole.pole_length*pole_dis_sin - cart_pole.pole_mass*cart_pole.pole_length*pole_dis_cos*(cart_pole.pole_mass*cart_pole.pole_length*pow(pole_vel,2.0)*pole_dis_sin - cart_pole.friction*cart_vel)) - cart_pole.pole_mass*cart_pole.pole_length*pole_dis_cos*(1.0/p)*force;

        cart_vel += (cart_accel + cart_accel_prev) / 2.0 * (time - time_prev);
        cart_dis += (cart_vel + cart_vel_prev) / 2.0 * (time - time_prev);

        pole_vel += (pole_accel + pole_accel_prev) / 2.0 * (time - time_prev);
        pole_dis += (pole_vel + pole_vel_prev) / 2.0 * (time - time_prev);

        outfile << time << "," << cart_accel << "," << cart_vel << "," << cart_dis << "," << pole_accel << "," << pole_vel << "," << pole_dis << "," <<  force << std::endl;

        /* == eigen decomposition method (Arimoto-Potter algorithm) == */
        clock_t start = clock();
        solveRiccatiArimotoPotter(A, B, Q, R, P);
        
        K_continous = R.inverse() * B.transpose()*P;
        // PRINT_MAT(K_continous);

        K_discrete = (R + B.transpose()*P*B).inverse()*B.transpose()*P*A;
        // PRINT_MAT(K_discrete);

        clock_t end = clock();
        // std::cout << "Computation time = " << (double)(end - start) / CLOCKS_PER_SEC << " (s)" << std::endl;

        X(0) = cart_dis;
        X(1) = cart_vel;
        X(2) = pole_dis;
        X(3) = pole_vel;

        u = -K_discrete * (X);
        force = u(0);

        if (time > 30.0){
            iteration = end_iteration;
        }

        cart_accel_prev = cart_accel;
        cart_vel_prev = cart_vel;
        cart_dis_prev = cart_dis;
        pole_accel_prev = pole_accel;
        pole_vel_prev = pole_vel;
        pole_dis_prev = pole_dis;

        time_prev = time;
        time = time + time_delta;
        iteration++;
    }
};
