#include "experiment.cpp"
#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <string>
#include "constraints.hpp"
#include <bitset>



// local_w anpassen 0.5 auf connected, kleiner auf 

using namespace std;
using namespace Eigen;


int main(){
    double band{2};
    ConstraintCoeffs<2> c1;
    c1.constype = "linear";
    c1.type = "eq";
    c1.coeffs = Vector2d{-1,-1};
    c1.cons = -1;
     //rrt_exp.run();

    std::vector<double> lb{-2.,-2.};
    std::vector<double> ub{2.,2.};

    Experiment<2,1> line("biased","uniform","","");
    line.setLocalBounds(lb, ub);
    line.setGlobalBounds(lb, ub);
    line.addConstraints(c1);
    line.setGlobalNumIter(5000);
    line.setLocalNumIter(0);
    line.setSave(true);
    line.evaluateKdesForPlot(false);
    line.run();

    return 1; 
} 
