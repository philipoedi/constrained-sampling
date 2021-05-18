#include <vector>
#include <iostream>
#include "optimizer.hpp"
#include "kde.hpp"
#include <string>
#include <Eigen/Dense>
#include "utils.hpp"
#include <cassert>

using namespace Eigen;

///////////////////////////////
// 1. Problem spcifications////
///////////////////////////////

std::string name_suffix = "diagonal_2";

// dimension - number of decision variables
const std::size_t n{2};
// number of inequalitiy constraints
const std::size_t m{1};
// number of equality constraints
const std::size_t l{0};

// method to for sampling
// {"biased", "slack"}
const std::string method{"slack"};

// lower bounds
const std::vector<double> lb{0,0};
// upper bounds
const std::vector<double> ub{1,1};

// linear constraint coefficients
// c.T@x - b

const double cs[m][n]{
    {-1,-1}};
const double bs[m]{-1};
std::vector<double> c{-1,-1};
double b{-1};



////////////////////
// 2. Simulations///
////////////////////

// number of samples
const int n_iter{1000};


////////////////////
// 3. Evaluation ///
////////////////////

// bandwidth estimator {"silverman", "scott"}
const std::string band_est{"silverman"};
// grid spacing
const double step{0.05};



int main()
{

   // 1. problem
    constraint_coeffs<n> ineqc;
    for (int i=0; i<n; i++)
    {
        ineqc.coeffs(i) = c[i];
    };
    ineqc.cons = b;
    ineqc.type = "ineq";
    ineqc.constype = "linear";

    kernel_estimator<n_iter,n> kdest;
    //std::vector<std::vector<double>> res;
    base_optimizer<n> opti();
    assert (method == "biased" || method == "slack");
    std::string name = utils::get_date_string()+"_"+ method + "_linear"+"_"+name_suffix;
    if (method == "biased")
    {
        biased_optimizer<n> opti(ineqc, lb, ub);
        opti.run(n_iter);
        //opti.results(res);
        opti.save_results(name+"_results");
        opti.save_samples(name+"_samples");
        kdest.fit(opti.results());
        kdest.find_optimal_bandwidth(band_est);
        kdest.predict(lb, ub, step);
        kdest.save_pdes(name+"_pdes");
    }
    else 
    {
        slack_optimizer<n,m,l> opti(ineqc, lb, ub);
        opti.run(n_iter);
        //opti.results(res);
        opti.save_results(name+"_results");
        opti.save_samples(name+"_samples");
        kdest.fit(opti.results());
        kdest.find_optimal_bandwidth(band_est);
        kdest.predict(lb, ub, step);
        kdest.save_pdes(name+"_pdes");
    };
     // saving results

    // writing resúlts to file
    //utils::write_vec2file(res,name+"_results");
    // writing samples/seeds to file
    //utils::write_vec2file(opti.samples(), name+"_seeds");
    // writing probability densities to file
    //utils::write_vec2file(,name);
    // writnig metadata
    utils::write_metadata2file("linear_cons.cpp",name);

   // ndvector2file(res, samples_name);
  /*opt.results(results);
    kdest.fit(ressults);
    kdest.predict(lb,ub,step);
*/
    return 0;
}
