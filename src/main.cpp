#include "sampler.hpp"
#include "optimizer.hpp"
#include <Eigen/Dense>
#include <iostream>
#include <nlopt.hpp>
#include <vector>


int main(int argc, char** argv){
    
    using namespace Eigen;
    using namespace nlopt;

    // dimension of dicision variable
    const std::size_t n = 2;
    const std::vector<double> lb{1., 2.};
    const std::vector<double> ub{5., 6.};
    Vector2d lb_v(lb.data());
    Vector2d ub_v(ub.data());

    // sampler setup
    UniformSampler<n> u(lb_v, ub_v);
    
    // optimizer setup
    Bias<n> b;
    ConstraintCoeffs<n> c;
    c.coeffs << 1., 1. ;
    c.cons = 6.;
    opt Biased("AUGLAG_EQ", n);
   // opt Biased("LD_MMA", n);
    opt local_opt("LD_SLSQP", n);
    local_opt.set_xtol_rel(1e-4);
    Biased.set_local_optimizer(local_opt);
    Biased.set_upper_bounds(ub);
    Biased.set_lower_bounds(lb);
    Biased.add_inequality_constraint(linearConstraint<n>, &c, 1e-8);
    Biased.set_xtol_rel(1e-4);

    // loop random sample optimize
    const int n_iter = 100000;
    std::vector<double> x(n);
    double minf;
    for(int i = 0; i < n_iter; ++i){
        b.x0 = u.sample();
  //      std::cout << "seed: " << b.x0[0] << ", " << b.x0[1] << std::endl;
        for (std::size_t i = 0; i < b.x0.size(); ++i) {
            x[i] = b.x0[i];
        }
        Biased.set_min_objective(Biased_objective<n>, &b);
        try{
            Biased.optimize(x, minf);
//            std::cout << "f(x): " << minf << std::endl;
            std::cout << minf << ", "  << x[0] << ", " << x[1] << ", "<<b.x0[0]  <<", "<< b.x0[1] << std::endl;
        }
        catch(std::exception &e){
            std::cerr << e.what() << std::endl;
        }
    }
    return 0;
}

