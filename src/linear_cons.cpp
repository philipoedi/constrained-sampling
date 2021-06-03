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

std::string name_suffix = "island";

// dimension - number of decision variables
const std::size_t n{2};
// number of inequalitiy constraints
const std::size_t m{2};
// number of equality constraints
const std::size_t l{0};

// method to for sampling
// {"Biased", "slack"}
const std::string method{"slack"};

// lower bounds
const std::vector<double> lb{0,0};
// upper bounds
const std::vector<double> ub{2,2};

// linear constraint coefficients
// c.T@x - b

const double cs[m][n] ={ 
    {-1,-1},
    {1,1}};
const double bs[m]{-1,2};
std::vector<double> c{-1,-1};
double b{-1};

std::vector<double> c2{1,1};
double b2{2};


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
    


   /* for (int i=0; i<m; i++)
    {
        std::cout << cs[i][0] << " "<< cs[i][1]  << std::endl;
    };*/

    ConstraintCoeffs<n> ineqc;
    for (int i=0; i<n; i++)
    {
        ineqc.coeffs(i) = c[i];
    };
    ineqc.cons = b;
    ineqc.type = "ineq";
    ineqc.constype = "linear";

    ConstraintCoeffs<n> ineqc2;
    for (int i=0; i<n; i++)
    {
        ineqc2.coeffs(i) = c2[i];
    };
    ineqc2.cons = b2;
    ineqc2.type = "ineq";
    ineqc2.constype = "linear";

    KernelEstimator<n_iter,n> kdest;
    //std::vector<std::vector<double>> res;
    BaseOptimizer<n> opti();
    assert (method == "Biased" || method == "slack");
    std::string name = utils::getDateString()+"_"+ method + "_linear"+"_"+name_suffix;
    if (method == "Biased")
    {
        BiasedOptimizer<n> opti(ineqc, lb, ub);
        opti.addConstraints(ineqc2);
        opti.run(n_iter);
        //opti.results(res);
        opti.saveResults(name+"_results");
        opti.save_samples(name+"_samples");
        kdest.fit(opti.results());
        kdest.find_optimal_bandwidth(band_est);
        kdest.predict(lb, ub, step);
        kdest.savePdes(name+"_pdes");
    }
    else 
    {
        SlackOptimizer<n,m,l> opti(ineqc2, lb, ub);
        opti.addConstraints(ineqc);
        opti.run(n_iter);
        //opti.results(res);
        opti.saveResults(name+"_results");
        opti.save_samples(name+"_samples");
        kdest.fit(opti.results());
        kdest.find_optimal_bandwidth(band_est);
        kdest.predict(lb, ub, step);
        kdest.savePdes(name+"_pdes");
    };
     // saving results

    // writing resúlts to file
    //utils::writeVec2File(res,name+"_results");
    // writing samples/seeds to file
    //utils::writeVec2File(opti.samples(), name+"_seeds");
    // writing probability densities to file
    //utils::writeVec2File(,name);
    // writnig metadata
    utils::writeMetadata2File("linear_cons.cpp",name);

   // ndvector2file(res, samples_name);
  /*opt.results(results);
    kdest.fit(ressults);
    kdest.predict(lb,ub,step);
*/
    return 0;
}
