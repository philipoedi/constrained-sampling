#include <iostream>
#include <string>
#include <Eigen/Dense>
#include "optimizer.hpp"
#include "sampler.hpp"
#include "constraints.hpp"


int main(){
    using namespace Eigen;
    const std::size_t n{2};

    std::vector<double> lb(n,0);
    std::vector<double> ub(n,2);
    std::cout << ub[0] << " " << ub[1] << std::endl;
    std::vector<double> w(n,0.5);

    ConstraintCoeffs<n> con;
    con.coeffs = Matrix<double,n,1>::Constant(-1);
    con.cons = -1 ;
    con.type = "ineq";
    con.constype = "linear"; 
    std::vector<ConstraintCoeffs<n>> cons;
    cons.push_back(con);

    GridWalk<n> gw;
    GridWalk<n> gw2(lb, ub, w);
    gw2.addConstraints(cons);
    std::cout << gw2.getBounds().first << std::endl;
    std::cout << gw2.getBounds().second << std::endl;
    gw2.run(5);
    gw2.saveResults("res");
    gw2.saveSamples("samples");

    std::cout << "ballwalk" << std::endl;
    gw2.makeBallWalk(0.5);
    gw2.run(5);
}
