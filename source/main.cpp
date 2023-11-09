#include <cstdlib>
#include <iomanip>
#include <iostream>

#include "simulation.hpp"

int main(int argc, char *argv[])
{
    std::cout << "argc == " << argc << '\n';
 
    for (int ndx{}; ndx != argc; ++ndx)
        std::cout << "argv[" << ndx << "] == " << std::quoted(argv[ndx]) << '\n';
    std::cout << "argv[" << argc << "] == "
              << static_cast<void*>(argv[argc]) << '\n';
 
    /*...*/

    Simulation simulation = Simulation();
    simulation.start(100);



    return argc == 3 ? EXIT_SUCCESS : EXIT_FAILURE; // optional return value
}