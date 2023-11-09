#ifndef CART_POLE_HPP
#define CART_POLE_HPP

class CartPole{
    public: 
        CartPole();
        ~CartPole();

        double cart_mass;

        double friction;
        double force;

        double pole_mass;
        double pole_length;
        double gravity = -9.81;
};
#endif