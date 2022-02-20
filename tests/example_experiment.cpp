#include "experiment.cpp"
#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <string>
#include "constraints.hpp"
#include <bitset>




using namespace std;
using namespace Eigen;


int main(){

    const int sphere_n{3};
    vector<double> local_ub_sphere_gw{1,1,1};
    vector<double> local_lb_sphere_gw{-1,-1,-1};
    vector<double> local_ub_sphere{0.5,0.5,0.5};
    vector<double> local_lb_sphere{-0.5,-0.5,-0.5};
    vector<double> local_lb_sphere_rrt{-0.25,-0.25,-0.25};
    vector<double> local_ub_sphere_rrt{0.25,0.25,0.25};
    vector<double> global_lb_sphere{-3,-2,-4};
    vector<double> global_ub_sphere{4,3,2};
    ConstraintCoeffs<sphere_n> sphere = createSphere<sphere_n>(1); // creates sphere-constraint with radius 1
    vector<double> local_w_2{0.2,0.2,0.2};
    vector<double> local_w{0.5,0.5,0.5};

    double filter_thresh = 0.4;
    //vector<double> local_w{0.5,0.5,0.5};

    // quadratic constraints that create disconnected manifold
    ConstraintCoeffs<3> y_z_cons = createConstraint2in3d<3>(1,2,-10,-1,-1.2,"ineq");
    ConstraintCoeffs<3> z_y_cons = createConstraint2in3d<3>(2,1,-10,-1,-1.2,"ineq");
    ConstraintCoeffs<3> x_z_cons = createConstraint2in3d<3>(0,2,-200,-1,-2,"ineq");


    // sphere center
    vector<double> global_lb_sphere_sym{-2,-2,-2};
    vector<double> global_ub_sphere_sym{2,2,2};
    // sphere off center
    vector<double> global_lb_sphere_non_sym{-3,-2,-4};
    vector<double> global_ub_sphere_non_sym{4,3,2};

    const double grid{0.025};
    const double band{5}; // bandwidth is band = 1/bw 
     
    // RRT
    // Experiment<n,m>, n-dim configuration space, m-eq-constraints
    Experiment<sphere_n,1> rrt("biased","uniform","biased","RRT");
    rrt.setGlobalBounds(global_lb_sphere, global_ub_sphere); // set the bounds of the configuration space
    rrt.setLocalBounds(local_lb_sphere_rrt, local_ub_sphere_rrt); // max size of the RRT 
    rrt.addConstraints(sphere); 
    // comment out below if closed manifold is needed
    rrt.addConstraints(y_z_cons);  
    rrt.addConstraints(z_y_cons); 
    rrt.addConstraints(x_z_cons); 
    // set 0 when 0 filter, or comment out
    rrt.setFilter(0.5);
    // set true when bandit sampler is applied
    rrt.setWeightedMultiple(true);
    // set true when on global stage rejection sampling appplied
    // if false then projection on to edges possible
    rrt.setGlobalEqOnly(false);
    rrt.setLocalEqOnly(false);
    // iterations/samples
    rrt.setGlobalNumIter(100);
    rrt.setLocalNumIter(50);
    rrt.setLocalAlpha(0.01); // alpha/stepsize of rrt
    rrt.setLocalUseTangent(true); // always set to true when equality constraints used
    rrt.setBandwidth(band);
    // for plotting of the probability density surface, KDE will be evaluated over a grid 
    rrt.setSphere(1); // specify 1 when equailty constraint is a sphere of radius 1, KDE will then be evaluated over the sphere surface not over full grid of configuration space
    rrt.setGridSpacing(0.5); // specifies how dense the grid is to evaluate the KDE for 	
    rrt.evaluateKdesForPlot(false);// can be set to false if many experiments are run and not all necessary to plot
    rrt.setSave(true);
   
    // gridwalk multiple 
    Experiment<sphere_n,1> gw_multi("biased","uniform","biased","grid-walk");
    gw_multi.setLocalBounds(local_lb_sphere_gw, local_ub_sphere_gw);
    gw_multi.addConstraints(sphere);
    gw_multi.setLocalNumIter(50);
    gw_multi.setGlobalNumIter(100);
    gw_multi.setLocalAlpha(0.005);
    gw_multi.setLocalWidths(local_w);
    gw_multi.setBandwidth(band);
    gw_multi.setSphere(1);
    gw_multi.setGridSpacing(0.5);
    gw_multi.setLocalUseTangent(true);
    gw_multi.setSave(true);
   
return 1;
}

