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
   
    ConstraintCoeffs<2> c1;
    c1.constype = "linear";
    c1.type = "ineq";
    c1.coeffs = Vector2d{-1,-1};
    c1.cons = -1;


    ConstraintCoeffs<1> x_cons_low = createConstraintParallelAxis<1>(0,1,0,"ineq");
    ConstraintCoeffs<1> x_cons_up = createConstraintParallelAxis<1>(0,-1,1,"ineq");


    vector<double> global_lb{-2};
    vector<double> global_ub{2};
    vector<double> local_lb{-1};
    vector<double> local_ub{1};
    vector<vector<double>> widths;
    widths.push_back(vector<double>{0.5});
    widths.push_back(vector<double>{0.4});
    widths.push_back(vector<double>{0.3});
    widths.push_back(vector<double>{0.2});
    widths.push_back(vector<double>{0.1});
    
    vector<vector<double>> seeds;
    vector<double> start{0.5};
    seeds.push_back(start);


    vector<vector<double>> ubs;
    vector<vector<double>> lbs;

    lbs.push_back(vector<double>{-1.1});
    ubs.push_back(vector<double>{1.1});
    lbs.push_back(vector<double>{-1.2});
    ubs.push_back(vector<double>{1.2});
    lbs.push_back(vector<double>{-1.3});
    ubs.push_back(vector<double>{1.3});
    lbs.push_back(vector<double>{-1.4});
    ubs.push_back(vector<double>{1.4});
    lbs.push_back(vector<double>{-1.5});
    ubs.push_back(vector<double>{1.5});


    const std::size_t n_samples{5000};
    const double grid{0.001};
    const double band{5};
    const std::size_t num_iter{100};
    Experiment<1,0> line_sampler("biased","uniform","biased","grid-walk");
    line_sampler.setLocalBounds(local_lb, local_ub);
    line_sampler.setGlobalBounds(global_lb, global_ub);
    line_sampler.addConstraints(x_cons_up);
    line_sampler.addConstraints(x_cons_low);
    line_sampler.setGridSpacing(0.025);
    line_sampler.setSave(true);
    line_sampler.setBandwidth(band);
    line_sampler.setGlobalEqOnly(true);
    line_sampler.setLocalEqOnly(true);
    line_sampler.setGlobalNumIter(1);
    line_sampler.setLocalNumIter(n_samples);
    line_sampler.setStartSeeds(seeds);
    
    Experiment<1,0> line_rrt("biased","uniform","biased","RRT");
    line_rrt.setLocalBounds(local_lb, local_ub);
    line_rrt.setGlobalBounds(global_lb, global_ub);
    line_rrt.addConstraints(x_cons_up);
    line_rrt.addConstraints(x_cons_low);
    line_rrt.setGridSpacing(0.025);
    line_rrt.setSave(true);
    line_rrt.setLocalAlpha(0.0025);
    line_rrt.setLocalNumIter(n_samples);
    line_rrt.setGlobalNumIter(1);
    line_rrt.setBandwidth(band);
    line_rrt.setLocalEqOnly(true);
    line_rrt.setStartSeeds(seeds);

    Experiment<1,0> ref("biased","uniform","","");
    ref.setLocalBounds(local_lb,local_ub);
    ref.setGlobalBounds(vector<double>{0},vector<double>{1});
    ref.setSave(true);
    ref.setGlobalNumIter(n_samples);
    ref.setLocalNumIter(0);
    ref.setBandwidth(band);
    rrt.saveKdesForPlot(false);

    for (std::size_t i=0; i < num_iter; i++){
        std::cout << widths.size() << std::endl;
        for (std::size_t j =0; j < widths.size(); j++){
            line_sampler.setLocalWidths(widths[j]);
            line_sampler.setLocalEqOnly(true);
            line_sampler.setSuffix("_eq_"+to_string(i)+"_w"+to_string(widths[j][0]));
            line_sampler.run();
            line_sampler.setLocalEqOnly(false);
            line_sampler.setSuffix("_ineq_"+to_string(i)+"_w"+to_string(widths[j][0]));
            line_sampler.run();
        }
        for (std::size_t l=0; l<ubs.size(); l++){
            line_rrt.setLocalBounds(lbs[l],ubs[l]);
            line_rrt.setLocalEqOnly(true);
            line_rrt.setSuffix("_eq_"+to_string(i)+"_b"+to_string(ubs[l][0]));
            line_rrt.run();
            line_rrt.setLocalEqOnly(false);
            line_rrt.setSuffix("_ineq_"+to_string(i)+"_b"+to_string(ubs[l][0]));
            line_rrt.run();
        }
        ref.setSuffix("_reference_"+to_string(i));
        ref.run();
    }

    return 1;
}
