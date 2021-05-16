#include "plotting.hpp"

void plot_line(Gnuplot &gp )
{
    std::cout << " plotting started" << std::endl;
    gp << "plot sin(x)\n"; 
}
