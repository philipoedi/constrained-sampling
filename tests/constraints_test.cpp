#include <iostream>
#include <cassert>
#include <Eigen/Dense>
#include <vector>
#include "constraints.hpp"
#include <functional>

using namespace Eigen;
using namespace std;
using namespace std::placeholders;
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
    auto f = std::bind(evaluateConstraint<n>, _1, c1);
    //std::function<double(std::vector<double>&)>;
    cout << "checkin bind\n" << f(x) << endl << endl;

    std::vector<double> grad = numericGradient(x, f, 1e-8);
    cout << "checking numeric gradient\n "<< grad[0] << " " << grad[1] << endl;
    std::vector<std::function<double(std::vector<double>&)>> funcs;
    funcs.push_back(f);
    funcs.push_back(f);
    Matrix<double,2,2> jac;
    jac = numericJacobian<2,2>(x,funcs,1e-8);
    cout << "numeric Jacobian\n" << jac << endl;

    for (int i=0;i<2;i++){
        funcs.push_back(std::bind(evaluateConstraint<n>,_1,c1));
    }
   
    Matrix<double,4,2> jac2;
    jac2 = numericJacobian<2,4>(x,funcs,1e-8);
    cout << "numeric Jacobian\n" << jac2 << endl;

    ConstraintCoeffs<3> sphere = createSphere<3>(1);
    std::vector<double> sphere_test{0,0,1};
    cout << "eval sphere at 0,0,1: " << evaluateConstraint<3>(sphere_test ,sphere ) << endl; 
    
    funcs.clear(),
    funcs.push_back(std::bind(evaluateConstraint<3>,_1,sphere));
    Matrix<double,1,3> jac3;
    jac3 = numericJacobian<3,1>(sphere_test, funcs, 1e-8);
    cout << "numeric Jacobian\n" << jac3 << endl;

};
