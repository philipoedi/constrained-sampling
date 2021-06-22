#include "experiment.cpp"
#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <string>
#include "constraints.hpp"

using namespace std;
using namespace Eigen;


int main(){
    const string global_optimizer{"biased"};
    const string global_sampler{"uniform"};    
    const string local_optimizer{"biased"};    
    const string local_sampler{"uniform"};    
    
    Experiment<2,2> rrt_exp;
    rrt_exp.setGlobalOptimizer(global_optimizer);
    rrt_exp.setLocalOptimizer(local_optimizer);
    rrt_exp.setGlobalSampler(global_sampler);
    rrt_exp.setLocalSampler(local_sampler);

    Experiment<2,2> rrt_exp2(
        global_optimizer,
        global_sampler,
        local_optimizer,
        local_sampler
    ); 
    
    vector<double> lb{0,0};
    vector<double> ub{2,2};
    rrt_exp.setLocalBounds(lb,ub);
    rrt_exp.setGlobalBounds(lb,ub);
   
    ConstraintCoeffs<2> c1;
    c1.constype = "linear";
    c1.type = "ineq";
    c1.coeffs = Vector2d{-1,-1};
    c1.cons = -1;
    cout << "coeffs" << c1.coeffs << endl;
    ConstraintCoeffs<2> c2;
    c2.constype = "linear";
    c2.type = "ineq";

    vector<ConstraintCoeffs<2>> cons;
    cons.push_back(c1);
    cons.push_back(c2);

    rrt_exp2.addConstraints(c1);
    rrt_exp2.addConstraints(cons);
    rrt_exp.addConstraints(c1);
    double local_alpha{2};
    double global_alpha{4};

    rrt_exp.setLocalAlpha(local_alpha);
    rrt_exp.setGlobalAlpha(global_alpha);
    
    rrt_exp.setGlobalNumIter(5);
    rrt_exp.setLocalNumIter(5);
    rrt_exp.run();

}
