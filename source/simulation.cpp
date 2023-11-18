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
    // std::ofstream outfile;
    // outfile.open(id + ".csv");
    this->id = id;

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



    CartPole cart_pole;

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

    // outfile <<  "time" << "," << "cart_vel_linear" << "," << "cart_dis_linear" <<  "," << "pole_vel_linear" << "," << "pole_dis_linear" << "," << "cart_accel" << "," << "cart_vel" << "," << "cart_dis" << "," << "pole_accel" << "," << "pole_vel" << "," << "pole_dis" << "," <<  "force" << std::endl;
    // outfile << time << "," <<  X(1) << "," << X(0) << "," << X(3) << "," << X(2) << "," << cart_accel << "," << cart_vel << "," <<  cart_dis << "," << pole_accel << "," << pole_vel << "," << pole_dis << "," << force << std::endl;
    // outfile.precision(12);

    double pole_dis_prev = pole_dis;
    double cart_dis_prev = cart_dis;
    double pole_vel_prev = pole_vel;
    double cart_vel_prev = cart_vel;
    double pole_accel_prev = pole_accel;
    double cart_accel_prev = cart_accel;

    double time_delta = 0.001;
    double time_prev = 0.0;


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


        this->times.push_back(time);
        this->cart_positions.push_back(cart_dis);
        this->cart_velocities.push_back(cart_vel);
        this->cart_accelerations.push_back(cart_accel);
        this->pole_positions.push_back(pole_dis);
        this->pole_velocities.push_back(pole_vel);
        this->pole_accelerations.push_back(pole_accel);
        this->forces.push_back(force);
        
        // outfile << time << "," <<  X(1) << "," << X(0) << "," << X(3) << "," << X(2) << "," << cart_accel << "," << cart_vel << "," <<  cart_dis << "," << pole_accel << "," << pole_vel << "," << pole_dis << "," << force << std::endl;

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
    std::map<std::string, std::vector<double>> key_value;
    this->simulation_data.data.try_emplace(this->id, key_value);
    this->simulation_data.data.at(this->id).try_emplace(std::string("time"), this->times);
    this->simulation_data.data.at(this->id).try_emplace(std::string("cart_displacement"), this->cart_positions);
    this->simulation_data.data.at(this->id).try_emplace(std::string("cart_velocity"), this->cart_velocities);
    this->simulation_data.data.at(this->id).try_emplace(std::string("cart_acceleration"), this->cart_accelerations);
    this->simulation_data.data.at(this->id).try_emplace(std::string("pole_displacement"), this->pole_positions);
    this->simulation_data.data.at(this->id).try_emplace(std::string("pole_velocity"), this->pole_velocities);
    this->simulation_data.data.at(this->id).try_emplace(std::string("pole_acceleration"), this->pole_accelerations);
    this->simulation_data.data.at(this->id).try_emplace(std::string("force"), this->forces);


    double decimal_point_accuracy = 0.05;
    if (X.isApprox(desired_state, decimal_point_accuracy)){
        this->status = 1.0;
        this->simulation_data.status.trajectory_status = 1.0;
    }
    this->simulation_data.status.trajectory_status = 0.0;

    // boolean will be used to catch errors one day...
    return true;
};

SimulationData Simulation::get_data(){
    return this->simulation_data;
}

namespace simulation_functions{
    // simple output function
    void output(std::string outfile, const std::vector<SimulationData>& m){
        std::ofstream out;
        out.open(outfile);
        out << "{" << std::endl;
        // Iterate using C++17 facilities
        for (int i_simulation = 0; i_simulation < m.size(); i_simulation++){
            int i_key_one = 0;
            for (const auto& [key_one, value_one] : m[i_simulation].data){
                if (i_key_one != m.size()){
                    out << "\"" << key_one << "\"" << ":" << "{" << std::endl;
                    out << "\"" << "status" << "\"" << ":" << m[i_simulation].status.trajectory_status <<  "," << std::endl;    
                    int i_key_two = 0;
                    for (const auto& [key_two, value_two] : value_one){
                        if (i_key_two != m.size()-1){
                            out << "\"" << key_two << "\"" << ":" << std::endl << "[";
                            for (int i_value = 0; i_value < value_two.size(); i_value++){
                                if (i_value != value_two.size()-1){
                                    out << value_two[i_value] << "," << " ";
                                }
                                else{
                                    out << value_two[i_value];
                                }
                            }
                            out << "]";
                            if (i_key_two != value_one.size()-1){
                                out << "," << std::endl;
                            }
                            else{

                            }
                        }
                        else{
                            out << "\"" << key_two << "\"" << ":" << std::endl;
                        }
                        i_key_two++;           
                    }
                    if (i_key_one != m.size()-1){
                        out << "}" << "," << std::endl;
                    }
                    else{
                        out << "}" << std::endl;
                    }
                }
                else{
                    out << "\"" << key_one << "\"" << ":" << "{";        
                }
                i_key_one++;
            }
        }
        out << '\n' << "}";
    }

}

