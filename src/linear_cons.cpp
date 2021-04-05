#include <vector>
#include <iostream>
#include "optimizer.hpp"
#include "kde.hpp"
#include <string>
#include <Eigen/Dense>

using namespace Eigen;

///////////////////////////////
// 1. Problem spcifications////
///////////////////////////////

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
const std::vector<double> ub{2,2};

// linear constraint coefficients
// c.T@x - b
std::vector<double> c{1,1};
double b{1};

////////////////////
// 2. Simulations///
////////////////////

// number of samples
const int n_iter{1000};


////////////////////
// 3. Evaluation ///
////////////////////

// bandwidth estimator {"silverman", "scott"}
const std::string band_est{"scott"};
// grid spacing
const double step{0.1};

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
    std::vector<std::vector<double>> res;

    if (method == "biased")
    {
        biased_optimizer<n> opti(ineqc, lb, ub);
        opti.run(n_iter);
        opti.results(res);
        kdest.fit(res);
        kdest.find_optimal_bandwidth(band_est);
        kdest.predict(lb, ub, step);

    }
    else if(method == "slack")
    {
        slack_optimizer<n,m,l> opti(ineqc, lb, ub);
        opti.run(n_iter);
        opti.results(res);
        kdest.fit(res);
        kdest.find_optimal_bandwidth(band_est);
        kdest.predict(lb, ub, step);
    }
    else
    {
        std::cout << "choose method 'slack' or 'biased'" << std::endl;
    };
       /*opt.results(results);
    kdest.fit(results);
    kdest.predict(lb,ub,step);
*/
    return 0;
}
