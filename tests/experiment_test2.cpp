	#include "experiment.cpp"
#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <string>
#include "constraints.hpp"
#include <bitset>



// local_w anpassen 0.5 auf connected, kleiner auf 

using namespace std;
using namespace Eigen;


int main(){
    const string global_optimizer{"biased"};
    const string global_sampler{"uniform"};    
    const string local_optimizer{"biased"};    
    const string local_sampler{"uniform"};    
    
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

    double local_alpha{2};
    double global_alpha{4};
   //rrt_exp.run();


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

    /*gw_multi.setSuffix("_off_center_disconnected_multiple_filter_ch_bandit");
    gw_multi.run();
 */






   /* 

        u_b_r_n.evaluateKdesForPlot(true);



    for (int i=0; i<num_experiments; i++){
       u_b_r_n.setSuffix("_disconnected_"+experiment_bit+"_testing_-1");
        u_b_r_n3.setSuffix("_"+experiment_bit+"_testing_-1");
        u_b_r_b_t.setSuffix("_disconnected_"+experiment_bit+"_testing_-1");
        u_b_r_b_t2.setSuffix("_"+experiment_bit+"_testing_-1");
        u_b_r_n.run();
        u_b_r_n3.run();
        u_b_r_b_t.run();
        u_b_r_b_t2.run();
        u_b_r_n.evaluateKdesForPlot(false);
        u_b_r_n3.evaluateKdesForPlot(false);
        u_b_r_b_t.evaluateKdesForPlot(false);
        u_b_r_b_t2.evaluateKdesForPlot(false);
        for (int i=0; i<n_samples_s; i++){
            u_b_r_n.setSuffix("_disconnected_"+experiment_bit+"_testing_"+to_string(i));
            u_b_r_n3.setSuffix("_"+experiment_bit+"_testing_"+to_string(i));
            u_b_r_b_t.setSuffix("_disconnected_"+experiment_bit+"_testing_"+to_string(i));
            u_b_r_b_t2.setSuffix("_"+experiment_bit+"_testing_"+to_string(i));
            u_b_r_n.run();
            u_b_r_n3.run();
            u_b_r_b_t.run();
            u_b_r_b_t2.run();
        } 
        u_b_r_n.evaluateKdesForPlot(true);
        u_b_r_n3.evaluateKdesForPlot(true);
        u_b_r_b_t.evaluateKdesForPlot(true);
        u_b_r_b_t2.evaluateKdesForPlot(true);
    
    }
    
    // biased optimizer only
    Experiment<sphere_n,1> u_b("biased","uniform","","");
    u_b.setLocalBounds(local_lb_sphere, local_ub_sphere);
    //u_b.setGlobalBounds(global_lb_sphere, global_ub_sphere);
    u_b.addConstraints(sphere);
    u_b.setGlobalNumIter(5000);
    u_b.setLocalNumIter(0);
    u_b.setLocalAlpha(0.005);
    u_b.setLocalWidths(local_w);
    u_b.setBandwidth(band);
    u_b.setSphere(1);
    u_b.setGridSpacing(0.5);
    u_b.setLocalUseTangent(true);
    u_b.setSave(true);
    
    Experiment<sphere_n,1> u_b2("biased","uniform","","");
    u_b2.setLocalBounds(local_lb_sphere, local_ub_sphere);
    //u_b2.setGlobalBounds(global_lb_sphere, global_ub_sphere);
    u_b2.addConstraints(z_y_cons);
    u_b2.addConstraints(x_z_cons);
    u_b2.addConstraints(y_z_cons);
    u_b2.addConstraints(sphere);
    u_b2.setGlobalNumIter(5000);
    u_b2.setLocalNumIter(0);
    u_b2.setLocalAlpha(0.005);
    u_b2.setLocalWidths(local_w);
    u_b2.setBandwidth(band);
    u_b2.setSphere(1);
    u_b2.setSave(true);
    u_b2.setGridSpacing(0.5);
    u_b2.setLocalUseTangent(true);
    

    const int num_factors_biased_only{2};
    int num_experiments_biased_only = pow(2,num_factors_biased_only);
    string experiment_bit_biased_only;
    for (int i=0; i<num_experiments_biased_only; i++){
        experiment_bit_biased_only = bitset<num_factors_biased_only>(i).to_string();
        if (experiment_bit_biased_only[0] == '1'){
            u_b.setGlobalEqOnly(true);
            u_b2.setGlobalEqOnly(true);
        } else {
            u_b.setGlobalEqOnly(false);
            u_b2.setGlobalEqOnly(false);
        }
        if (experiment_bit_biased_only[1] == '1'){
            u_b.setGlobalBounds(global_lb_sphere_sym, global_ub_sphere_sym);
            u_b2.setGlobalBounds(global_lb_sphere_sym, global_ub_sphere_sym);
        } else {
            u_b.setGlobalBounds(global_lb_sphere_non_sym, global_ub_sphere_non_sym);
            u_b2.setGlobalBounds(global_lb_sphere_non_sym, global_ub_sphere_non_sym);
        }
        u_b2.setSuffix("_disconnected_"+experiment_bit_biased_only+"_testing_-1");
        u_b.setSuffix("_"+experiment_bit_biased_only+"_testing_-1");
        u_b.run();
        u_b2.run();
        u_b.evaluateKdesForPlot(false);
        u_b2.evaluateKdesForPlot(false);
        for (int i=0; i<n_samples_s; i++){
            u_b.setSuffix("_"+experiment_bit_biased_only+"_testing_"+to_string(i));
            u_b2.setSuffix("_disconnected_"+experiment_bit_biased_only+"_testing_"+to_string(i));
            u_b.run();
            u_b2.run();
        } 
        u_b.evaluateKdesForPlot(true);
        u_b2.evaluateKdesForPlot(true);
    
    }
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
    sphere_ref.setSuffix("_testing_-1");
    sphere_ref.run();


    Experiment<3,1> sphere_ref2("sphere","reference","","");
    sphere_ref2.setLocalBounds(local_lb_sphere, local_lb_sphere);
    sphere_ref2.setGlobalBounds(global_lb_sphere, global_ub_sphere);
    sphere_ref2.setLocalNumIter(0);
    sphere_ref2.setGlobalNumIter(5000);
    sphere_ref2.setSphere(1);
    sphere_ref2.setSave(true);
    sphere_ref2.setLocalUseTangent(false);
    sphere_ref2.setBandwidth(band);
    sphere_ref2.addConstraints(z_y_cons);
    sphere_ref2.addConstraints(x_z_cons);
    sphere_ref2.addConstraints(y_z_cons);
    sphere_ref2.setGridSpacing(0.5);
    sphere_ref2.setSuffix("_disconnnected_testing_-1");
    sphere_ref2.run();

    sphere_ref.evaluateKdesForPlot(false);
    sphere_ref2.evaluateKdesForPlot(false);

    for (int i=0; i<n_samples_s; i++){
        sphere_ref.setSuffix("_testing_"+to_string(i));
        sphere_ref2.setSuffix("_disconnected_testing_"+to_string(i));
        sphere_ref.run();
        sphere_ref2.run();
    }


    return 1;
    
    for (int i=0; i<n_samples_s; i++){
    //sphere_ref.setSuffix("_disconnected_ref_testing_"+to_string(i));
        //u_b_r_n.setSuffix("_disconnected_filter_testing_"+to_string(i));
        u_b.setSuffix("_disconnected_testing_"+to_string(i));
        u_b_r_n2.setSuffix("_disconnected_no_filter_testing_"+to_string(i));
        u_b_r_n3.setSuffix("_disconnected_no_fitler_wide_testing_"+to_string(i));
        u_b_r_n4.setSuffix("_disconnected_filter_wide_testing_"+to_string(i));
        u_b_r_b_t.setSuffix("_disconnected_multiple_testing_"+to_string(i));
        u_b_r_b_t2.setSuffix("_disconnected_single_testing_"+to_string(i));
        sphere_ref.setSuffix("_disconnected_testing_"+to_string(i));
        u_b_r_n.run();
        u_b.run();
        u_b_r_n2.run();
        u_b_r_n3.run();
        u_b_r_n4.run();
        u_b_r_b_t.run();
        u_b_r_b_t2.run();
        sphere_ref.run();
    } 


    return 1;

    
    // gridwalk
   

   // gridwalk single global
    Experiment<sphere_n,1> u_b_r_b_t3("biased","uniform","biased","grid-walk");
    u_b_r_b_t3.setLocalBounds(local_lb_sphere, local_ub_sphere);
    u_b_r_b_t3.setGlobalBounds(global_lb_sphere, global_ub_sphere);
    u_b_r_b_t3.addConstraints(sphere);
    u_b_r_b_t3.setGlobalNumIter(1);
    u_b_r_b_t3.setLocalNumIter(5000);
    u_b_r_b_t3.setLocalAlpha(0.005);
    u_b_r_b_t3.setLocalWidths(local_w);
    u_b_r_b_t3.setBandwidth(band);
    u_b_r_b_t3.setSphere(1);
    u_b_r_b_t3.setGridSpacing(0.5);
    u_b_r_b_t3.setLocalUseTangent(true);
    u_b_r_b_t3.setSave(true);
    u_b_r_b_t3.setSuffix("_single");
    u_b_r_b_t3.run();
    
    // RRT many global no filter
    Experiment<sphere_n,1> u_b_r_n2("biased","uniform","biased","RRT");
    u_b_r_n2.setLocalBounds(local_lb_sphere, local_ub_sphere);
    u_b_r_n2.setGlobalBounds(global_lb_sphere,global_ub_sphere);
    u_b_r_n2.addConstraints(sphere);
    u_b_r_n2.setGlobalNumIter(100);
    u_b_r_n2.setLocalNumIter(50);
    u_b_r_n2.setLocalAlpha(0.005);
    u_b_r_n2.setLocalUseTangent(true);
    u_b_r_n2.setBandwidth(band);
    u_b_r_n2.setSphere(1);
    u_b_r_n2.setGridSpacing(0.5);
    u_b_r_n2.setSave(true);
    u_b_r_n2.setSuffix("_no_filter");
    u_b_r_n2.run();
 
    
    return 1;
    

    // RRT fe global, wide local filter 
    Experiment<sphere_n,1> u_b_r_n4("biased","uniform","biased","RRT");
    u_b_r_n4.setLocalBounds(local_lb_sphere3, local_ub_sphere3);
    u_b_r_n4.setGlobalBounds(global_lb_sphere,global_ub_sphere);
    u_b_r_n4.addConstraints(sphere);
    u_b_r_n4.setGlobalNumIter(100);
    u_b_r_n4.setLocalNumIter(50);
    u_b_r_n4.setLocalAlpha(0.005);
    u_b_r_n4.setLocalUseTangent(true);
    u_b_r_n4.setBandwidth(band);
    u_b_r_n4.setSphere(1);
    u_b_r_n4.setFilter(0.25);
    u_b_r_n4.setGridSpacing(0.5);
    u_b_r_n4.setSave(true);
    u_b_r_n4.setSuffix("_filter_wide");
    u_b_r_n4.run();



   
    */

    // 2d examples

    /*
    const int n2d{2};
    //std::vector<double> local_bounds_lb_2d{-0.5,-0.5};
    //std::vector<double> local_bounds_ub_2d{0.5,0.5};
    //std::vector<double> global_bounds_lb_2d{0,0};
    //std::vector<double> global_bounds_ub_2d{2,2};


    Experiment<n2d,1> u_b_2d("biased","uniform","","");
    u_b_2d.setLocalBounds(local_lb_circle, local_ub_circle);
    u_b_2d.setGlobalBounds(global_lb_circle, global_ub_circle);
    u_b_2d.addConstraints(c1);
    u_b_2d.setLocalNumIter(0);
    u_b_2d.setGlobalNumIter(5000);
    u_b_2d.setBandwidth(band);
    u_b_2d.setGridSpacing(grid);
    u_b_2d.setLocalUseTangent(false);
    u_b_2d.setSave(true);
    u_b_2d.setSuffix("_2d");
    u_b_2d.run();
    */
/*
    // reference sphere
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
    sphere_ref.run();
 */   // reference circle
    /*
    Experiment<2,1> circle_ref("circle","reference","","");
    circle_ref.setLocalBounds(local_lb_circle, local_ub_circle);
    circle_ref.setGlobalBounds(global_lb_circle, global_ub_circle);
    circle_ref.setSphere(1);
    circle_ref.setLocalNumIter(0);
    circle_ref.setGlobalNumIter(5000);
    circle_ref.setSave(true);
    circle_ref.setBandwidth(band);
    circle_ref.setGridSpacing(grid);
    //circle_ref.run();

    // reference line
    Experiment<2,1> line_ref("line","reference","","");
    line_ref.setLocalBounds(local_lb_circle, local_ub_circle);
    line_ref.setGlobalBounds(global_lb_circle, global_ub_circle);
    line_ref.addConstraints(c1); 
    line_ref.setLocalNumIter(0);
    line_ref.setGlobalNumIter(5000);
    line_ref.setSave(true);
    line_ref.setBandwidth(band);
    line_ref.setGridSpacing(grid);
    //line_ref.run();

    // reference ineq linear
    Experiment<2,1> ineq_ref("rejection","reference","","");
    ineq_ref.setLocalBounds(local_lb_circle, local_ub_circle);
    ineq_ref.setGlobalBounds(global_lb_circle,global_ub_circle);
    ineq_ref.addConstraints(c1);
    ineq_ref.setLocalNumIter(0);
    ineq_ref.setGlobalNumIter(100);
    ineq_ref.setSave(true);
    ineq_ref.setBandwidth(band);
    ineq_ref.setGridSpacing(grid);
    //ineq_ref.run();


    std::vector<int> num_samples{5000};
    std::vector<double> num_band{2};

    std::cout << "----------------------" << std::endl;
    std::cout << "Referencens" << std::endl;
    std::cout << "----------------------" << std::endl;
    int current_num_samples;
    double current_band;
    std::string suf;
    for (int i=0; i<num_samples.size(); i++){
        for (int j=0; j<num_samples.size(); j++){
            current_num_samples = num_samples[j];
            current_band = num_band[i];
            suf = std::to_string(current_num_samples) + "-" + std::to_string(current_band);
            ineq_ref.setGlobalNumIter(current_num_samples);
            line_ref.setGlobalNumIter(current_num_samples);
            circle_ref.setGlobalNumIter(current_num_samples);
            sphere_ref.setGlobalNumIter(current_num_samples);
            ineq_ref.setBandwidth(current_band);
            line_ref.setBandwidth(current_band);
            circle_ref.setBandwidth(current_band);
            sphere_ref.setBandwidth(current_band);
            ineq_ref.setSuffix(suf);
            line_ref.setSuffix(suf);
            circle_ref.setSuffix(suf);
            sphere_ref.setSuffix(suf);
            ineq_ref.run();
            line_ref.run();
            circle_ref.run();
            sphere_ref.run();
        }
    }
    
    std::cout << "----------------------" << std::endl;
    std::cout << "2d gridwalk" << std::endl;
    std::cout << "----------------------" << std::endl;
     
    // gridwalk
    Experiment<2,1> u_b_2d_gw("biased","uniform","","grid-walk");
    u_b_2d_gw.setLocalBounds(local_lb_circle, local_ub_circle);
    u_b_2d_gw.setGlobalBounds(global_lb_circle, global_ub_circle);
    u_b_2d_gw.addConstraints(c1);
    u_b_2d_gw.setLocalWidths(local_ub_2d);
    u_b_2d_gw.setLocalNumIter(5000);
    u_b_2d_gw.setGlobalNumIter(1);
    u_b_2d_gw.setBandwidth(band);
    u_b_2d_gw.setGridSpacing(grid);
    u_b_2d_gw.setLocalUseTangent(false);
    u_b_2d_gw.setSave(true);
    u_b_2d_gw.setSuffix("_2d");
    u_b_2d_gw.run();
    
    std::cout << "----------------------" << std::endl;
    std::cout << "2d rrt" << std::endl;
    std::cout << "----------------------" << std::endl;
     
    // gridwalk
    Experiment<2,1> u_b_2d_rrt("biased","uniform","","RRT");
    u_b_2d_rrt.setLocalBounds(local_lb_circle, local_ub_circle);
    u_b_2d_rrt.setGlobalBounds(global_lb_circle, global_ub_circle);
    u_b_2d_rrt.addConstraints(c1);
    u_b_2d_rrt.setLocalWidths(local_ub_2d);
    u_b_2d_rrt.setLocalAlpha(0.05);
    u_b_2d_rrt.setLocalNumIter(5000);
    u_b_2d_rrt.setGlobalNumIter(1);
    u_b_2d_rrt.setBandwidth(band);
    u_b_2d_rrt.setGridSpacing(grid);
    u_b_2d_rrt.setLocalUseTangent(false);
    u_b_2d_rrt.setSave(true);
    u_b_2d_rrt.setSuffix("_2d");
    u_b_2d_rrt.run();


    std::cout << "----------------------" << std::endl;
    std::cout << "biased equality" << std::endl;
    std::cout << "----------------------" << std::endl;
 

    // based uniform only eq 
    c1.type = "eq";
    Experiment<2,1> u_b_2d_eq("biased","uniform","","");
    u_b_2d_eq.setLocalBounds(local_lb_circle, local_ub_circle);
    u_b_2d_eq.setGlobalBounds(global_lb_circle, global_ub_circle);
    u_b_2d_eq.addConstraints(c1);
    u_b_2d_eq.setLocalNumIter(0);
    u_b_2d_eq.setGlobalNumIter(5000);
    u_b_2d_eq.setBandwidth(band);
    u_b_2d_eq.setGridSpacing(grid);
    u_b_2d_eq.setLocalUseTangent(false);
    u_b_2d_eq.setSave(true);
    u_b_2d_eq.setSuffix("_2d_eq");
    u_b_2d_eq.run();

    // biased uniform circle
    std::cout << "----------------------" << std::endl;
    std::cout << "circle equality" << std::endl;
    std::cout << "----------------------" << std::endl;
 

    ConstraintCoeffs<2> circle = createSphere<2>(1);
    Experiment<2,1> u_b_2d_circle("biased","uniform","","");
    u_b_2d_circle.setLocalBounds(local_lb_circle, local_ub_circle);
    u_b_2d_circle.setGlobalBounds(global_lb_circle, global_ub_circle);
    u_b_2d_circle.addConstraints(circle);
    u_b_2d_circle.setLocalNumIter(0);
    u_b_2d_circle.setGlobalNumIter(5000);
    u_b_2d_circle.setBandwidth(band);
    u_b_2d_circle.setGridSpacing(grid);
    u_b_2d_circle.setLocalUseTangent(false);
    u_b_2d_circle.setSave(true);
    u_b_2d_circle.setSuffix("_2d_circle");
    u_b_2d_circle.run();

*/

    // ------------------
    // Multiple samples
    // ------------------

    /*
    for (int i=0; i<n_samples_s; i++){
        sphere_ref.setSuffix("_ref_testing_"+to_string(i));
        u_b_r_n.setSuffix("_filter_testing_"+to_string(i));
        u_b_r_n2.setSuffix("_no_filter_testing_"+to_string(i));
        u_b_r_n3.setSuffix("_no_fitler_wide_testing_"+to_string(i));
        u_b_r_n4.setSuffix("_filter_wide_testing_"+to_string(i));
        u_b_r_b_t.setSuffix("_multiple_testing_"+to_string(i));
        u_b_r_b_t2.setSuffix("_single_testing_"+to_string(i));
        u_b.setSuffix("_testing_"+to_string(i));
        sphere_ref.setSuffix("_testing_"+to_string(i));
        u_b_r_n.run();
        u_b_r_n2.run();
        u_b_r_n3.run();
        u_b_r_n4.run();
        u_b_r_b_t.run();
        u_b_r_b_t2.run();
        u_b.run();
        sphere_ref.run();
    } 
*/ /*
    sphere_ref.addConstraints(y_z_cons);
    sphere_ref.addConstraints(z_y_cons);
    sphere_ref.addConstraints(x_z_cons);
    sphere_ref.setSuffix("_disconnected");
    sphere_ref.evaluateKdesForPlot(true);
    sphere_ref.run();

    u_b_r_n.evaluateKdesForPlot(true);
    u_b_r_n.addConstraints(y_z_cons);
    u_b_r_n.addConstraints(z_y_cons);
    u_b_r_n.addConstraints(x_z_cons);
    u_b_r_n.setSuffix("_filter_disconnected");
    u_b_r_n.run();
    
    u_b_r_n2.evaluateKdesForPlot(true);
    u_b_r_n2.addConstraints(y_z_cons);
    u_b_r_n2.addConstraints(z_y_cons);
    u_b_r_n2.addConstraints(x_z_cons);
    u_b_r_n2.setSuffix("_no_filter_disconnected");
    u_b_r_n2.run();
    
    u_b_r_n3.evaluateKdesForPlot(true);
    u_b_r_n3.addConstraints(y_z_cons);
    u_b_r_n3.addConstraints(z_y_cons);
    u_b_r_n3.addConstraints(x_z_cons);
    u_b_r_n3.setSuffix("_no_filter_wide_disconnected");
    u_b_r_n3.run();
      
    u_b_r_n4.evaluateKdesForPlot(true);
    u_b_r_n4.addConstraints(y_z_cons);
    u_b_r_n4.addConstraints(z_y_cons);
    u_b_r_n4.addConstraints(x_z_cons);
    u_b_r_n4.setSuffix("_filter_wide_disconnected");
    u_b_r_n4.run();
     
    u_b_r_b_t.evaluateKdesForPlot(true);
    u_b_r_b_t.addConstraints(y_z_cons);
    u_b_r_b_t.addConstraints(z_y_cons);
    u_b_r_b_t.addConstraints(x_z_cons);
    u_b_r_b_t.setSuffix("_multiple_disconnected");
    u_b_r_b_t.run();
     
    u_b_r_b_t2.evaluateKdesForPlot(true);
    u_b_r_b_t2.addConstraints(y_z_cons);
    u_b_r_b_t2.addConstraints(z_y_cons);
    u_b_r_b_t2.addConstraints(x_z_cons);
    u_b_r_b_t2.setSuffix("_single_disconnected");
    u_b_r_b_t2.run();
   

    u_b.evaluateKdesForPlot(true);
    u_b.addConstraints(y_z_cons);
    u_b.addConstraints(z_y_cons);
    u_b.addConstraints(x_z_cons);
    u_b.setSuffix("_disconnected");
    u_b.run(); 


    sphere_ref.evaluateKdesForPlot(false);
    u_b_r_n.evaluateKdesForPlot(false);
    u_b_r_n2.evaluateKdesForPlot(false);
    u_b_r_n3.evaluateKdesForPlot(false);
    u_b_r_n4.evaluateKdesForPlot(false);
    u_b_r_b_t.evaluateKdesForPlot(false);
    u_b_r_b_t2.evaluateKdesForPlot(false);
    u_b.evaluateKdesForPlot(false);*/
/*
    for (int i=0; i<n_samples_s; i++){
        sphere_ref.setSuffix("_disconnected_ref_testing_"+to_string(i));
        u_b_r_n.setSuffix("_disconnected_filter_testing_"+to_string(i));
        u_b.setSuffix("_disconnected_testing_"+to_string(i));
        u_b_r_n2.setSuffix("_disconnected_no_filter_testing_"+to_string(i));
        u_b_r_n3.setSuffix("_disconnected_no_fitler_wide_testing_"+to_string(i));
        u_b_r_n4.setSuffix("_disconnected_filter_wide_testing_"+to_string(i));
        u_b_r_b_t.setSuffix("_disconnected_multiple_testing_"+to_string(i));
        u_b_r_b_t2.setSuffix("_disconnected_single_testing_"+to_string(i));
        sphere_ref.setSuffix("_disconnected_testing_"+to_string(i));
        u_b_r_n.run();
        u_b.run();
        u_b_r_n2.run();
        u_b_r_n3.run();
        u_b_r_n4.run();
        u_b_r_b_t.run();
        u_b_r_b_t2.run();
        sphere_ref.run();
    } 
*/

    // uniform_biased + metropolis_hastings_rejection

    // uniform_rejection + metropolis_hastings_biased

    // uniform_rejection + rrt_rejection
}
