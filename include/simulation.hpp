class Simulation{
    public:
        Simulation();
        ~Simulation();
        void start(unsigned int end_iteration);

    private:
        unsigned int iteration;
        unsigned int end_iteration;


};