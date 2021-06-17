#include "optimizer.hpp"
#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include "constraints.hpp"
#include <string>


using namespace Eigen;

int main(){
    const std::size_t n{2};

    ConstraintCoeffs<n> c1;
    ConstraintCoeffs<n> c2;
    c1.coeffs << 2, 0; 
    c1.cons = 0;
    c1.type = "eq";
    c1.constype = "linear";
    c2.constype= "quadratic";
    c2.type = "eq";
    c2.r = 0.5;
    c2.P = Matrix2d::Identity();
    std::vector<double> lb{-2,-2};
    std::vector<double> ub{+2,+2};
    BiasedOptimizer<n> bopt;
    bopt.setBounds(lb,ub);
    bopt.addConstraints(c1);
    bopt.addConstraints(c2);
    const int n_iter{1};
    bopt.run(n_iter);
    std::vector<std::vector<double>> results;
    bopt.results(results);
    bopt.saveResults("tangent");
};
