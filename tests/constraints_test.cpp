#include <iostream>
#include <cassert>
#include <Eigen/Dense>
#include <vector>
#include "constraints.hpp"

using namespace Eigen;
using namespace std;

int main(){
    const std::size_t n{2};
    ConstraintCoeffs<n> c1;
    c1.constype = "linear";
    c1.type = "ineq";
    c1.coeffs = Vector2d{-1,-1};
    c1.cons = -1;

    ConstraintCoeffs<n> c2;
    c2.constype = "linear";
    c2.type = "ineq";
    c2.coeffs = Vector2d{-1,1};
    c2.cons = 1;

    vector<double> x{0,0};
    cout << "Evaluate constraint: " << endl;
    cout << evaluateConstraint<n>(x, c1) << endl;
    
    cout << "Is feasible: " << endl;
    cout << isFeasible<n>(x,c1) << endl;

    std::vector<ConstraintCoeffs<n>> cons;
    cons.push_back(c1);
    cons.push_back(c2);

    cout << "Check multiple constraints feasible: " << endl;
    cout << isFeasibleM<n>(x, cons) << endl;
    
    Vector2d xe{1,1};
    Vector2d lb{0,0};
    Vector2d ub{2,2};
    cout << "Check within bounds"<< endl;
    cout << boundsCheck<n>(xe,lb,ub) <<endl;
};
