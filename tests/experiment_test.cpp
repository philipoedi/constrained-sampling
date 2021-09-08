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

    double local_alpha{2};
    double global_alpha{4};
   //rrt_exp.run();


    const int sphere_n{3};
    vector<double> local_ub_sphere{0.5,0.5,0.5};
    vector<double> local_lb_sphere{-0.5,-0.5,-0.5};
    vector<double> global_lb_sphere{-3,-2,-4};
    vector<double> global_ub_sphere{4,3,2};
    ConstraintCoeffs<sphere_n> sphere = createSphere<sphere_n>(1);
    vector<double> local_w{0.5,0.5,0.5};

    //constraint on x axis
    ConstraintCoeffs<3> x_axis_cons = createConstraintParallelAxis<3>(0,1,0.2,"ineq");
    //ConstraintCoeffs<3> x_axis_cons2 = createConstraintParallelAxis<3>(0,1,-0.2,"ineq");
    ConstraintCoeffs<3> y_z_cons = createConstraint2in3d<3>(1,2,-10,-1,-1.2,"ineq");
    ConstraintCoeffs<3> z_y_cons = createConstraint2in3d<3>(2,1,-10,-1,-1.2,"ineq");
    ConstraintCoeffs<3> x_z_cons = createConstraint2in3d<3>(0,2,-200,-1,-2,"ineq");


    vector<double> global_lb_circle{-2,-2};
    vector<double> global_ub_circle{2,2};
    vector<double> local_lb_circle{-2,-2};
    vector<double> local_ub_circle{-2,-2};
    vector<double> local_lb_2d{-0.5,-0.5};
    vector<double> local_ub_2d{0.5,0.5};

    const double grid{0.025};
    const double band{2};
    const int local_n_iter{100};
    const int global_n_iter{200};
    
    // RRT many global and filter 
    Experiment<sphere_n,1> u_b_r_n("biased","uniform","biased","RRT");
    u_b_r_n.setLocalBounds(local_lb_sphere, local_ub_sphere);
    u_b_r_n.setGlobalBounds(global_lb_sphere,global_ub_sphere);
    u_b_r_n.addConstraints(sphere);
    //u_b_r_n.addConstraints(x_axis_cons);
    u_b_r_n.addConstraints(y_z_cons);
    u_b_r_n.addConstraints(z_y_cons);
    u_b_r_n.addConstraints(x_z_cons);
    u_b_r_n.setGlobalNumIter(100);
    u_b_r_n.setLocalNumIter(50);
    u_b_r_n.setLocalAlpha(0.01);
    u_b_r_n.setLocalUseTangent(true);
    u_b_r_n.setBandwidth(band);
    u_b_r_n.setSphere(1);
    u_b_r_n.setGridSpacing(0.5);
    u_b_r_n.setFilter(0.25);
    u_b_r_n.setSave(true);
    u_b_r_n.setSuffix("_filter_XXXXXX");
    u_b_r_n.run();  

    return 1;
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
    //u_b_r_n2.run();
 
    
    // RRT fe global, wide local 
    vector<double> local_lb_sphere3{-1,-1,-1};
    vector<double> local_ub_sphere3{1,1,1};
    Experiment<sphere_n,1> u_b_r_n3("biased","uniform","biased","RRT");
    u_b_r_n3.setLocalBounds(local_lb_sphere3, local_ub_sphere3);
    u_b_r_n3.setGlobalBounds(global_lb_sphere,global_ub_sphere);
    u_b_r_n3.addConstraints(sphere);
    u_b_r_n3.setGlobalNumIter(100);
    u_b_r_n3.setLocalNumIter(50);
    u_b_r_n3.setLocalAlpha(0.005);
    u_b_r_n3.setLocalUseTangent(true);
    u_b_r_n3.setBandwidth(band);
    u_b_r_n3.setSphere(1);
    u_b_r_n3.setGridSpacing(0.5);
    u_b_r_n3.setSave(true);
    u_b_r_n3.setSuffix("_no_filter_wide");
    //u_b_r_n3.run();

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
    //u_b_r_n4.run();



    // gridwalk
    Experiment<sphere_n,1> u_b_r_b_t("biased","uniform","biased","grid-walk");
    u_b_r_b_t.setLocalBounds(local_lb_sphere, local_ub_sphere);
    u_b_r_b_t.setGlobalBounds(global_lb_sphere, global_ub_sphere);
    u_b_r_b_t.addConstraints(sphere);
    u_b_r_b_t.setGlobalNumIter(100);
    u_b_r_b_t.setLocalNumIter(50);
    u_b_r_b_t.setLocalAlpha(0.005);
    u_b_r_b_t.setLocalWidths(local_w);
    u_b_r_b_t.setBandwidth(band);
    u_b_r_b_t.setSphere(1);
    u_b_r_b_t.setGridSpacing(0.5);
    u_b_r_b_t.setLocalUseTangent(true);
    u_b_r_b_t.setSave(true);
    u_b_r_b_t.setSuffix("_multiple");
    u_b_r_b_t.run();
 

   // gridwalk single global
    Experiment<sphere_n,1> u_b_r_b_t2("biased","uniform","biased","grid-walk");
    u_b_r_b_t2.setLocalBounds(local_lb_sphere, local_ub_sphere);
    u_b_r_b_t2.setGlobalBounds(global_lb_sphere, global_ub_sphere);
    u_b_r_b_t2.addConstraints(sphere);
    u_b_r_b_t2.setGlobalNumIter(1);
    u_b_r_b_t2.setLocalNumIter(5000);
    u_b_r_b_t2.setLocalAlpha(0.005);
    u_b_r_b_t2.setLocalWidths(local_w);
    u_b_r_b_t2.setBandwidth(band);
    u_b_r_b_t2.setSphere(1);
    u_b_r_b_t2.setGridSpacing(0.5);
    u_b_r_b_t2.setLocalUseTangent(true);
    u_b_r_b_t2.setSave(true);
    u_b_r_b_t2.setSuffix("_single");
    u_b_r_b_t2.run();
    
    // biased optimizer only
    Experiment<sphere_n,1> u_b("biased","uniform","","");
    u_b.setLocalBounds(local_lb_sphere, local_ub_sphere);
    u_b.setGlobalBounds(global_lb_sphere, global_ub_sphere);
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
    u_b.setSuffix("");
    u_b.run();
    


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
    // reference circle
    
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

    const int n_samples_s{100};


    sphere_ref.evaluateKdesForPlot(false);
    u_b_r_n.evaluateKdesForPlot(false);
    u_b_r_n2.evaluateKdesForPlot(false);
    u_b_r_n3.evaluateKdesForPlot(false);
    u_b_r_n4.evaluateKdesForPlot(false);
    u_b_r_b_t.evaluateKdesForPlot(false);
    u_b_r_b_t2.evaluateKdesForPlot(false);
    u_b.evaluateKdesForPlot(false);

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

    // uniform_biased + metropolis_hastings_rejection

    // uniform_rejection + metropolis_hastings_biased

    // uniform_rejection + rrt_rejection



}
