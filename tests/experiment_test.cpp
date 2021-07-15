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
    
    Experiment<2,1> rrt_exp;
    rrt_exp.setGlobalOptimizer(global_optimizer);
    rrt_exp.setLocalOptimizer(local_optimizer);
    rrt_exp.setGlobalSampler(global_sampler);
    rrt_exp.setLocalSampler(local_sampler);

    Experiment<2,1> rrt_exp2(
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
    //rrt_exp.run();
 
    const int sphere_n{3};
    Experiment<sphere_n,1> u_b_r_n("biased","uniform","biased","RRT");
    vector<double> local_lb_sphere{-0.25,-0.25,-0.25};
    vector<double> local_ub_sphere{0.25,0.25,0.25};
    u_b_r_n.setLocalBounds(local_lb_sphere, local_ub_sphere);
    vector<double> global_lb_sphere{-2,-2,-2};
    vector<double> global_ub_sphere{2,2,2};
    u_b_r_n.setGlobalBounds(global_lb_sphere,global_ub_sphere);
    ConstraintCoeffs<sphere_n> sphere = createSphere<sphere_n>(1);
    u_b_r_n.addConstraints(sphere);
    u_b_r_n.setGlobalNumIter(5);
    u_b_r_n.setLocalNumIter(10);
    u_b_r_n.setLocalAlpha(0.005);
    u_b_r_n.setLocalUseTangent(true);
    u_b_r_n.setBandwidth(5);
    u_b_r_n.setSphere(1);
    u_b_r_n.setGridSpacing(0.5);
    u_b_r_n.setSave(true);
    u_b_r_n.run();



    Experiment<2,1> u_b_r_b_t("biased","uniform","biased","grid-walk");
    vector<double> local_lb{-0.25,-0.25};
    vector<double> local_ub{0.25,0.25};
    u_b_r_b_t.setLocalBounds(local_lb, local_ub);
    u_b_r_n.setLocalBounds(local_lb, local_ub);
    u_b_r_b_t.setGlobalBounds(lb,ub);
    u_b_r_b_t.addConstraints(c1);
    u_b_r_b_t.setGlobalNumIter(5);
    u_b_r_b_t.setLocalNumIter(5);
    u_b_r_b_t.setLocalAlpha(0.005);
    u_b_r_b_t.run();
    

    // uniform_biased + metropolis_hastings_rejection

    // uniform_rejection + metropolis_hastings_biased

    // uniform_rejection + rrt_rejection



}
