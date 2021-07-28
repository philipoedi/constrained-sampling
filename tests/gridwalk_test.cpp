#include <iostream>
#include <string>
#include <Eigen/Dense>
#include "optimizer.hpp"
#include "sampler.hpp"
#include "constraints.hpp"


int main(){
    using namespace Eigen;
    const std::size_t n{2};
    const std::size_t m{1};
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

    ConstraintCoeffs<n> con2;
    con.coeffs = Matrix<double,n,1>::Constant(-1);
    con.cons = -1 ;
    con.type = "eq";
    con.constype = "linear"; 
    std::vector<ConstraintCoeffs<n>> cons2;

    GridWalk<n,m> gw;
    GridWalk<n,m> gw2(lb, ub, w);
    gw2.addConstraints(cons);
    std::cout << gw2.getBounds().first << std::endl;
    std::cout << gw2.getBounds().second << std::endl;
    gw2.run(5);
    gw2.saveResults("res");
    std::cout << "results saved" << std::endl;
    gw2.saveSamples("samples");
    std::cout << "samples saved" << std::endl;
    std::cout << "ballwalk" << std::endl;
    gw2.makeBallWalk(0.5);
    gw2.run(5);
    GridWalk<n,m> gw3(lb, ub,w);
    gw3.addConstraints(cons2);
    std::cout << "has optimizer: " << gw2.hasOptimizer() << std::endl;
    gw3.setOptimizer("biased", lb, ub);
    std::cout << "has optimizer: " <<  gw2.hasOptimizer() << std::endl;
    gw3.run(5);
    gw3.setRunOnTangent(true);
    std::cout << "running on tangent" << std::endl;
    gw3.runOnTangent(5);
    gw3.setRunOnTangent(false);
    gw3.runOnTangent(5);
    std::cout << "checking constructors " << std::endl;
    /*
    GridWalk<n,m> gw3(lb, ub, w, true);
    GridWalk<n,m> gw4(lb, ub, w, true, 0.3);
    gw3.addConstraints(cons);
    gw3.setOptimizer("biased", lb, ub);
    gw3.run(5);
    gw4.addConstraints(cons);
    gw4.setOptimizer("biased", lb, ub);
    
    gw4.run(5);
    */
}
