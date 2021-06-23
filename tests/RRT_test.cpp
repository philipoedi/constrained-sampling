#include <vector>
#include <Eigen/Dense>
#include "RRT.hpp"


int main(){
    using namespace Eigen;
    const std::size_t n{2};
    std::vector<double> vec1{1,2};
    Vector2d eig1(vec1.data()); 
    node<n> n1;
    n1.location = eig1;
    std::cout << n1() << std::endl; 
    
    // bounds
    std::vector<double> lb{0,0};
    std::vector<double> ub{2,2};
    
    // step size
    double alpha{0.5};

    // create RRT obj
    RRT<n> rrt;
    rrt.setBounds(lb, ub);
    rrt.setStepSize(alpha);

    // add constraints
    ConstraintCoeffs<n> c;
    c.type = "ineq";
    c.constype = "linear";
    std::vector<double> coeff_vec{-1,-1};
    Matrix<double,n,1> coeff(coeff_vec.data());
    c.coeffs = coeff;
    c.cons = -1;         
    std::vector<ConstraintCoeffs<n>> cons;
    cons.push_back(c);
    rrt.addConstraints(cons);
    // run
    const int n_iter{100};
    rrt.run(n_iter);
    // save files
    rrt.saveResults("results");
    rrt.saveSamples("samples"); 

    std::vector<double> seed{1,1};
    rrt.runOnTangent(3, seed, lb, ub); 

    return 0;
};
