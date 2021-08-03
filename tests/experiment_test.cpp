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
    vector<double> local_ub_sphere{0.5,0.5,0.5};
    vector<double> local_lb_sphere{-0.5,-0.5,-0.5};
    vector<double> global_lb_sphere{-2,-2,-2};
    vector<double> global_ub_sphere{2,2,2};
    ConstraintCoeffs<sphere_n> sphere = createSphere<sphere_n>(1);
    vector<double> local_w{0.5,0.5,0.5};

    const double band{5};
    const int local_yn_iter{100};
    const int global_n_iter{200};
/*
    // RRT many global and filter 
    Experiment<sphere_n,1> u_b_r_n("biased","uniform","biased","RRT");
    u_b_r_n.setLocalBounds(local_lb_sphere, local_ub_sphere);
    u_b_r_n.setGlobalBounds(global_lb_sphere,global_ub_sphere);
    u_b_r_n.addConstraints(sphere);
    u_b_r_n.setGlobalNumIter(200);
    u_b_r_n.setLocalNumIter(100);
    u_b_r_n.setLocalAlpha(0.01);
    u_b_r_n.setLocalUseTangent(true);
    u_b_r_n.setBandwidth(band);
    u_b_r_n.setSphere(1);
    u_b_r_n.setGridSpacing(0.5);
    u_b_r_n.setFilter(0.25);
    u_b_r_n.setSave(true);
    u_b_r_n.setSuffix("_filter");
    u_b_r_n.run();


    // RRT many global no filter
    
    Experiment<sphere_n,1> u_b_r_n2("biased","uniform","biased","RRT");
    u_b_r_n2.setLocalBounds(local_lb_sphere, local_ub_sphere);
    u_b_r_n2.setGlobalBounds(global_lb_sphere,global_ub_sphere);
    u_b_r_n2.addConstraints(sphere);
    u_b_r_n2.setGlobalNumIter(200);
    u_b_r_n2.setLocalNumIter(100);
    u_b_r_n2.setLocalAlpha(0.005);
    u_b_r_n2.setLocalUseTangent(true);
    u_b_r_n2.setBandwidth(band);
    u_b_r_n2.setSphere(1);
    u_b_r_n2.setGridSpacing(0.5);
    u_b_r_n2.setSave(true);
    u_b_r_n2.setSuffix("_no_filter");
    u_b_r_n2.run();
 
    
    // RRT fe global, wide local 
    vector<double> local_lb_sphere3{-1,-1,-1};
    vector<double> local_ub_sphere3{1,1,1};
    Experiment<sphere_n,1> u_b_r_n3("biased","uniform","biased","RRT");
    u_b_r_n3.setLocalBounds(local_lb_sphere3, local_ub_sphere3);
    u_b_r_n3.setGlobalBounds(global_lb_sphere,global_ub_sphere);
    u_b_r_n3.addConstraints(sphere);
    u_b_r_n3.setGlobalNumIter(200);
    u_b_r_n3.setLocalNumIter(100);
    u_b_r_n3.setLocalAlpha(0.005);
    u_b_r_n3.setLocalUseTangent(true);
    u_b_r_n3.setBandwidth(band);
    u_b_r_n3.setSphere(1);
    u_b_r_n3.setGridSpacing(0.5);
    u_b_r_n3.setSave(true);
    u_b_r_n3.setSuffix("_no_filter_wide");
    u_b_r_n3.run();



    // gridwalk
    Experiment<sphere_n,1> u_b_r_b_t("biased","uniform","biased","grid-walk");
    u_b_r_b_t.setLocalBounds(local_lb_sphere, local_ub_sphere);
    u_b_r_b_t.setGlobalBounds(global_lb_sphere, global_ub_sphere);
    u_b_r_b_t.addConstraints(sphere);
    u_b_r_b_t.setGlobalNumIter(200);
    u_b_r_b_t.setLocalNumIter(100);
    u_b_r_b_t.setLocalAlpha(0.005);
    u_b_r_b_t.setLocalWidths(local_w);
    u_b_r_b_t.setBandwidth(band);
    u_b_r_b_t.setSphere(1);
    u_b_r_b_t.setGridSpacing(0.5);
    u_b_r_b_t.setLocalUseTangent(true);
    u_b_r_b_t.setSave(true);
    u_b_r_b_t.setSuffix("_multiple");
    u_b_r_b_t.run();
 */ 
   
/*
    // gridwalk single global
    Experiment<sphere_n,1> u_b_r_b_t2("biased","uniform","biased","grid-walk");
    u_b_r_b_t2.setLocalBounds(local_lb_sphere, local_ub_sphere);
    u_b_r_b_t2.setGlobalBounds(global_lb_sphere, global_ub_sphere);
    u_b_r_b_t2.addConstraints(sphere);
    u_b_r_b_t2.setGlobalNumIter(1);
    u_b_r_b_t2.setLocalNumIter(20000);
    u_b_r_b_t2.setLocalAlpha(0.005);
    u_b_r_b_t2.setLocalWidths(local_w);
    u_b_r_b_t2.setBandwidth(band);
    u_b_r_b_t2.setSphere(1);
    u_b_r_b_t2.setGridSpacing(0.5);
    u_b_r_b_t2.setLocalUseTangent(true);
    u_b_r_b_t2.setSave(true);
    u_b_r_b_t2.setSuffix("_single");
    u_b_r_b_t2.run();
*/

    // 2d examples

    const int n2d{2};
    std::vector<double> local_bounds_lb_2d{-0.5,-0.5};
    std::vector<double> local_bounds_ub_2d{0.5,0.5};
    std::vector<double> global_bounds_lb_2d{0,0};
    std::vector<double> global_bounds_ub_2d{2,2};


    Experiment<n2d,1> u_b_2d("biased","uniform","","");
    u_b_2d.setLocalBounds(local_bounds_lb_2d, local_bounds_ub_2d);
    u_b_2d.setGlobalBounds(global_bounds_lb_2d, global_bounds_ub_2d);
    u_b_2d.addConstraints(c1);
    u_b_2d.setLocalNumIter(0);
    u_b_2d.setGlobalNumIter(100);
    u_b_2d.setBandwidth(band);
    u_b_2d.setGridSpacing(0.1);
    u_b_2d.setLocalUseTangent(false);
    u_b_2d.setSave(true);
    u_b_2d.setSuffix("_2d");
    u_b_2d.run();


    

    // uniform_biased + metropolis_hastings_rejection

    // uniform_rejection + metropolis_hastings_biased

    // uniform_rejection + rrt_rejection



}
