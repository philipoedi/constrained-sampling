#include <vector>
#include <iostream>
#include "gnuplot-iostream.h"
#include "plotting.hpp"
#include <string>

int main()
{
    Gnuplot gp;
    //plot_line(gp);
    
    std::vector<std::vector<double>> data(5, std::vector<double>(3,3));
    gp << "splot '-' \n";
    gp.send2d(data);
}

