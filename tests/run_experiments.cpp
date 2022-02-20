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
    ConstraintCoeffs<sphere_n> sphere = createSphere<sphere_n>(1);
    vector<double> local_w_2{0.2,0.2,0.2};
    vector<double> local_w{0.5,0.5,0.5};

    double filter_thresh = 0.4;
    //vector<double> local_w{0.5,0.5,0.5};

    //constraint on x axis
    ConstraintCoeffs<3> x_axis_cons = createConstraintParallelAxis<3>(0,1,0.2,"ineq");
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
    const double band{5};
    const int local_n_iter{100};
    const int global_n_iter{200};
    
    // RRT
    Experiment<sphere_n,1> rrt("biased","uniform","biased","RRT");
    rrt.setLocalBounds(local_lb_sphere_rrt, local_ub_sphere_rrt);
    rrt.addConstraints(sphere);
    rrt.setGlobalNumIter(100);
    rrt.setLocalNumIter(50);
    rrt.setLocalAlpha(0.01);
    rrt.setLocalUseTangent(true);
    rrt.setBandwidth(band);
    rrt.setSphere(1);
    rrt.setGridSpacing(0.5);
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
   
    // gridwalk single 
    Experiment<sphere_n,1> gw_single("biased","uniform","biased","grid-walk");
    gw_single.setLocalBounds(local_lb_sphere_gw, local_ub_sphere_gw);
    gw_single.addConstraints(sphere);
    gw_single.setLocalNumIter(5000);
    gw_single.setGlobalNumIter(1);
    gw_single.setLocalAlpha(0.005);
    gw_single.setLocalWidths(local_w);
    gw_single.setBandwidth(band);
    gw_single.setSphere(1);
    gw_single.setGridSpacing(0.5);
    gw_single.setLocalUseTangent(true);
    gw_single.setSave(true);
   
    // biased optimizer only
    Experiment<sphere_n,1> biased("biased","uniform","","");
    biased.setLocalBounds(local_lb_sphere, local_ub_sphere);
    biased.addConstraints(sphere);
    biased.setGlobalNumIter(5000);
    biased.setLocalNumIter(0);
    biased.setLocalAlpha(0.005);
    biased.setLocalWidths(local_w);
    biased.setBandwidth(band);
    biased.setSphere(1);
    biased.setGridSpacing(0.5);
    biased.setLocalUseTangent(true);
    biased.setSave(true);
   
    // reference
    Experiment<3,1> sphere_ref("sphere","reference","","");
    sphere_ref.setLocalBounds(local_lb_sphere, local_lb_sphere);
    sphere_ref.setGlobalBounds(global_lb_sphere, global_ub_sphere);
    sphere_ref.setLocalNumIter(0);
    sphere_ref.setGlobalNumIter(5000);
    sphere_ref.setSphere(1);
    sphere_ref.setSave(true);
    sphere_ref.setLocalUseTangent(false);
    sphere_ref.setBandwidth(band);
    sphere_ref.setGridSpacing(0.5);

    const int n_samples{100};
           
    // center connected for plot
    biased.setSuffix("_center_connected_-1");
    gw_single.setSuffix("_center_connected_single_-1");
    gw_multi.setSuffix("_center_connected_multiple_-1");
    rrt.setSuffix("_center_connected_-1");
    sphere_ref.setSuffix("_center_connected_-1");
    biased.setGlobalBounds(global_lb_sphere_sym, global_ub_sphere_sym);
    gw_single.setGlobalBounds(global_lb_sphere_sym, global_ub_sphere_sym);
    gw_multi.setGlobalBounds(global_lb_sphere_sym, global_ub_sphere_sym);
    rrt.setGlobalBounds(global_lb_sphere_sym, global_ub_sphere_sym);
    biased.run();
    gw_single.run();
    //gw_multi.run();
    rrt.run();
    sphere_ref.run();
    // center connected for eval
    gw_single.evaluateKdesForPlot(false);
    biased.evaluateKdesForPlot(false);
    //gw_multi.evaluateKdesForPlot(false);
    rrt.evaluateKdesForPlot(false);
    sphere_ref.evaluateKdesForPlot(false);
   
    for (int i=0; i<n_samples; i++){
        biased.setSuffix("_center_connected_"+to_string(i));
        gw_single.setSuffix("_center_connected_single_"+to_string(i));
        //gw_multi.setSuffix("_center_connected_multiple_"+to_string(i));
        rrt.setSuffix("_center_connected_"+to_string(i));
        sphere_ref.setSuffix("_center_connected_"+to_string(i));
        biased.run();
        gw_single.run();
        //gw_multi.run();
        rrt.run();
        sphere_ref.run();
    }
    gw_single.evaluateKdesForPlot(true);
    biased.evaluateKdesForPlot(true);
    gw_multi.evaluateKdesForPlot(true);
    //rrt.evaluateKdesForPlot(true);
    sphere_ref.evaluateKdesForPlot(true);
   

    // off center connected
    biased.setSuffix("_off_center_connected_-1");
    gw_single.setSuffix("_off_center_connected_single_-1");
    gw_multi.setSuffix("_off_center_connected_multiple_-1");
    rrt.setSuffix("_off_center_connected_-1");
    sphere_ref.setSuffix("_off_center_connected_-1");
    biased.setGlobalBounds(global_lb_sphere_non_sym, global_ub_sphere_non_sym);
    //gw_single.setGlobalBounds(global_lb_sphere_non_sym, global_ub_sphere_non_sym);
    gw_multi.setGlobalBounds(global_lb_sphere_non_sym, global_ub_sphere_non_sym);
    rrt.setGlobalBounds(global_lb_sphere_non_sym, global_ub_sphere_non_sym);
    biased.run();
    //gw_single.run();
    gw_multi.run();
    rrt.run();
    sphere_ref.run();
 
    gw_multi.setSuffix("_off_center_connected_multiple_filter_-1");
    rrt.setSuffix("_off_center_connected_filter_-1");
    gw_multi.setFilter(filter_thresh);
    rrt.setFilter(filter_thresh);
    gw_multi.run();
    rrt.run();

    // off center connected for eval
    //gw_single.evaluateKdesForPlot(false);
    gw_multi.evaluateKdesForPlot(false);
    rrt.evaluateKdesForPlot(false);
    sphere_ref.evaluateKdesForPlot(false);
    biased.evaluateKdesForPlot(false);
   
    for (int i=0; i<n_samples; i++){
        gw_multi.setFilter(0);
        rrt.setFilter(0);
        biased.setSuffix("_off_center_connected_"+to_string(i));
        //gw_single.setSuffix("_off_center_connected_single_"+to_string(i));
        gw_multi.setSuffix("_off_center_connected_multiple_"+to_string(i));
        rrt.setSuffix("_off_center_connected_"+to_string(i));
        sphere_ref.setSuffix("_off_center_connected_"+to_string(i));
        biased.run();
        //gw_single.run();
        gw_multi.run();
        rrt.run();
        sphere_ref.run();

        gw_multi.setSuffix("_off_center_connected_multiple_filter_"+to_string(i));
        rrt.setSuffix("_off_center_connected_filter_"+to_string(i));
        gw_multi.setFilter(filter_thresh);
        rrt.setFilter(filter_thresh);
        gw_multi.run();
        rrt.run();
    }


    //gw_single.evaluateKdesForPlot(true);
    gw_multi.evaluateKdesForPlot(true);
    rrt.evaluateKdesForPlot(true);
    sphere_ref.evaluateKdesForPlot(true);
    biased.evaluateKdesForPlot(true);
    
   
    // off-center disconnected 
    rrt.setFilter(0);
    gw_multi.setFilter(0);

    biased.addConstraints(y_z_cons);
    biased.addConstraints(x_z_cons);
    biased.addConstraints(z_y_cons);
    rrt.addConstraints(y_z_cons);
    rrt.addConstraints(x_z_cons);
    rrt.addConstraints(z_y_cons);
    gw_multi.addConstraints(y_z_cons);
    gw_multi.addConstraints(x_z_cons);
    gw_multi.addConstraints(z_y_cons);
    gw_multi.setLocalWidths(local_w);
    //gw_single.addConstraints(y_z_cons);
    //gw_single.addConstraints(x_z_cons);
    //gw_single.addConstraints(z_y_cons);
    sphere_ref.addConstraints(y_z_cons);
    sphere_ref.addConstraints(x_z_cons);
    sphere_ref.addConstraints(z_y_cons);
    
    //biased.setSuffix("_off_center_disconnected_-1");
    //gw_single.setSuffix("_off_center_disconnected_single_-1");
    gw_multi.setSuffix("_off_center_disconnected_multiple_-1");
    rrt.setSuffix("_off_center_disconnected_-1");
    sphere_ref.setSuffix("_off_center_disconnected_-1");
    //biased.run();
    //gw_single.run();
    gw_multi.run();
    rrt.run();
    sphere_ref.run();


    gw_multi.setFilter(filter_thresh); 
    rrt.setFilter(filter_thresh); 
    gw_multi.setSuffix("_off_center_disconnected_multiple_filter_-1");
    rrt.setSuffix("_off_center_disconnected_filter_-1");
    gw_multi.run();
    rrt.run();


    gw_multi.setGlobalEqOnly(true);
    gw_multi.setLocalEqOnly(true);
    rrt.setGlobalEqOnly(true); 
    rrt.setLocalEqOnly(true); 
    gw_multi.setSuffix("_off_center_disconnected_multiple_filter_ch_-1");
    rrt.setSuffix("_off_center_disconnected_filter_ch_-1");
    gw_multi.run();
    rrt.run();

    gw_multi.setWeightedMultiple(true);
    gw_multi.setSuffix("_off_center_disconnected_multiple_filter_ch_bandit_-1");
    gw_multi.run();
    rrt.setWeightedMultiple(true);
    rrt.setSuffix("_off_center_disconnected_filter_ch_bandit_-1");
    rrt.run();
    //gw_single.evaluateKdesForPlot(false);
    gw_multi.evaluateKdesForPlot(false);
    rrt.evaluateKdesForPlot(false);
    sphere_ref.evaluateKdesForPlot(false);
   
    for (int i=0; i<n_samples; i++){
        gw_multi.setWeightedMultiple(false);
        gw_multi.setGlobalEqOnly(false);
        gw_multi.setLocalEqOnly(false);
        rrt.setWeightedMultiple(false);
        rrt.setLocalEqOnly(false);
        rrt.setGlobalEqOnly(false);
        gw_multi.setFilter(0);
        rrt.setFilter(0);
        //biased.setSuffix("_off_center_disconnected_"+to_string(i));
        //gw_single.setSuffix("_off_center_disconnected_single_"+to_string(i));
        gw_multi.setSuffix("_off_center_disconnected_multiple_"+to_string(i));
        rrt.setSuffix("_off_center_disconnected_"+to_string(i));
        sphere_ref.setSuffix("_off_center_disconnected_"+to_string(i));
        //biased.run();
        //gw_single.run();
        gw_multi.run();
        rrt.run();
        sphere_ref.run();

        gw_multi.setSuffix("_off_center_disconnected_multiple_filter_"+to_string(i));
        rrt.setSuffix("_off_center_disconnected_filter_"+to_string(i));
        gw_multi.setFilter(filter_thresh);
        rrt.setFilter(filter_thresh);
        gw_multi.run();
        rrt.run();
 
        gw_multi.setSuffix("_off_center_disconnected_multiple_filter_ch_"+to_string(i));
        rrt.setSuffix("_off_center_disconnected_filter_ch_"+to_string(i));
        gw_multi.setGlobalEqOnly(true);
        gw_multi.setLocalEqOnly(true);
        rrt.setLocalEqOnly(true);
        rrt.setGlobalEqOnly(true);
        gw_multi.run();
        rrt.run();

        gw_multi.setWeightedMultiple(true);
        gw_multi.setSuffix("_off_center_disconnected_multiple_filter_ch_bandit_"+to_string(i));
        gw_multi.run();
        rrt.setSuffix("_off_center_disconnected_filter_ch_bandit_"+to_string(i));
        rrt.setWeightedMultiple(true);
        rrt.run();
      
    }
return 1;
}

