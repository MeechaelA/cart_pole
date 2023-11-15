#include "simulation.hpp"

#define PRINT_MAT(X) std::cout << #X << ":\n" << X << std::endl << std::endl

Simulation::Simulation(){
    this->iteration = 0;
    this->end_iteration = 0;
};

Simulation::~Simulation(){};


void solveRiccati(Eigen::MatrixXd &A, Eigen::MatrixXd &B,
                  Eigen::MatrixXd &Q, Eigen::MatrixXd &R,
                  Eigen::MatrixXd &P, Eigen::MatrixXd &K, double dt)
{
    double tolerance = 0.00001;
    int max_iter = 100000;
    
    P = Q;
    Eigen::MatrixXd P_next;
    double diff = std::numeric_limits<double>::max();
    int ii = 0;

    while (diff > tolerance || ii > max_iter)
    {

        ii++;
        P_next = P + (A.transpose() * P + P * A + Q - P * B * R.inverse() * B.transpose() * P) * dt;
        diff = fabs((P_next - P).maxCoeff());
        P = P_next;
        // std::cout << diff << "\n";
    }

    K = R.inverse() * B.transpose() * P;
}

bool Simulation::start(std::string id, double end_time, unsigned int end_iteration, Eigen::MatrixXd Q, Eigen::MatrixXd R, Eigen::MatrixXd start_state, Eigen::MatrixXd desired_state){
    std::ofstream outfile;
    outfile.open(id + ".csv");

    this->end_iteration = end_iteration;

    uint dim_x = 4;
    uint dim_u = 1;

    Eigen::MatrixXd X = start_state;
    Eigen::MatrixXd X_prev = Eigen::MatrixXd::Zero(dim_x, dim_u);
    Eigen::MatrixXd X_dot = Eigen::MatrixXd::Zero(dim_x, dim_u);
    Eigen::MatrixXd X_dot_prev = Eigen::MatrixXd::Zero(dim_x, dim_u);
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(dim_x, dim_x);
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(dim_x, dim_u);
    Eigen::MatrixXd P = Eigen::MatrixXd::Zero(dim_x, dim_x);
    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(dim_x, dim_x);
    Eigen::MatrixXd K_discrete = Eigen::MatrixXd::Zero(dim_x, dim_x);
    Eigen::MatrixXd u = Eigen::MatrixXd::Zero(dim_u, dim_u);

    X_prev = X;

    X_dot(0) = 0.0;
    X_dot(1) = 0.0;
    X_dot(2) = 0.0;
    X_dot(3) = 0.0;
    X_dot_prev = X_dot;



    CartPole    cart_pole;

    cart_pole.friction = 1.0;
    cart_pole.cart_mass = 5.0;
    cart_pole.pole_mass = 1.0;
    cart_pole.pole_length = 1.0;
    double b = 1.0;
    
    double cart_dis = X(0);
    double cart_vel = X(1);
    double pole_dis = X(2);
    double pole_vel = X(3);
    double pole_accel = 0.0;
    double cart_accel = 0.0;

    double force = 0.0;
    double time = 0.0;

    outfile <<  "time" << "," << "cart_vel_linear" << "," << "cart_dis_linear" <<  "," << "pole_vel_linear" << "," << "pole_dis_linear" << "," << "cart_accel" << "," << "cart_vel" << "," << "cart_dis" << "," << "pole_accel" << "," << "pole_vel" << "," << "pole_dis" << "," <<  "force" << std::endl;
    outfile << time << "," <<  X(1) << "," << X(0) << "," << X(3) << "," << X(2) << "," << cart_accel << "," << cart_vel << "," <<  cart_dis << "," << pole_accel << "," << pole_vel << "," << pole_dis << "," << force << std::endl;

    double pole_dis_prev = pole_dis;
    double cart_dis_prev = cart_dis;
    double pole_vel_prev = pole_vel;
    double cart_vel_prev = cart_vel;
    double pole_accel_prev = pole_accel;
    double cart_accel_prev = cart_accel;

    double time_delta = 0.001;
    double time_prev = 0.0;

    outfile.precision(12);

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

    // clock_t start = clock();
    solveRiccati(A, B, Q, R, P, K, time_delta);
    // clock_t end = clock();

    // PRINT_MAT(K);

    while(this->iteration < this->end_iteration){
        double pole_dis_cos = cos(pole_dis);
        double pole_dis_sin = sin(pole_dis);
        // std::cout << "cos: "<< pole_dis_cos << std::endl;
        // std::cout << "sin: "<< pole_dis_cos << std::endl;
        double p = cart_pole.pole_mass * cart_pole.pole_length * cart_pole.pole_length * (cart_pole.cart_mass + cart_pole.pole_mass*pole_dis_cos);

        cart_accel = (1.0/p)*(pow(-cart_pole.pole_mass,2.0)*pow(cart_pole.pole_length,2.0)*cart_pole.gravity*pole_dis_cos*pole_dis_sin + cart_pole.pole_mass*pow(cart_pole.pole_length,2.0)*(cart_pole.pole_mass*cart_pole.pole_length*pow(pole_vel,2.0)*pole_dis_sin - cart_pole.friction*cart_vel)) + cart_pole.pole_mass*cart_pole.pole_length*cart_pole.pole_length*(1.0/p)*force;
        pole_accel = (1.0/p)*((cart_pole.pole_mass+cart_pole.cart_mass)*cart_pole.pole_mass*cart_pole.gravity*cart_pole.pole_length*pole_dis_sin - cart_pole.pole_mass*cart_pole.pole_length*pole_dis_cos*(cart_pole.pole_mass*cart_pole.pole_length*pow(pole_vel,2.0)*pole_dis_sin - cart_pole.friction*cart_vel)) - cart_pole.pole_mass*cart_pole.pole_length*pole_dis_cos*(1.0/p)*force;

        cart_vel = (cart_accel + cart_accel_prev) / 2.0 * (time - time_prev);
        cart_dis = (cart_vel + cart_vel_prev) / 2.0 * (time - time_prev);

        pole_vel = ((pole_accel + pole_accel_prev) / 2.0 * (time - time_prev));
        pole_dis = ((pole_vel + pole_vel_prev) / 2.0 * (time - time_prev));

        u = -1.0*K * (X - desired_state);
        // X(0) = cart_dis;
        // X(1) = cart_vel;
        // X(2) = pole_dis;
        // X(3) = pole_vel;
        X = X + A * X * time_delta + B * u * time_delta;

        force = u(0);

        outfile << time << "," <<  X(1) << "," << X(0) << "," << X(3) << "," << X(2) << "," << cart_accel << "," << cart_vel << "," <<  cart_dis << "," << pole_accel << "," << pole_vel << "," << pole_dis << "," << force << std::endl;

        cart_accel_prev = cart_accel;
        cart_vel_prev = cart_vel;
        cart_dis_prev = cart_dis;
        pole_accel_prev = pole_accel;
        pole_vel_prev = pole_vel;
        pole_dis_prev = pole_dis;

        time_prev = time;
        time = time + time_delta;

        if (time > end_time){
            this->iteration = this->end_iteration;
        }
        this->iteration++;
    }

    double decimal_point_accuracy = 0.05;
    if (X.isApprox(desired_state, decimal_point_accuracy)){
        return true;
    }
    return false;
};
