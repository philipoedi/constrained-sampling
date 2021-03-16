#include "metrics.hpp"
#include <vector>
#include <iostream>
#include <cassert>
#include <Eigen/Dense>
#include <cmath>
#include <limits>
#include <nlopt.hpp>

using namespace Eigen;
using namespace nlopt;

int main(){
    double eps = std::numeric_limits<double>::epsilon();
    // test objective function coverage optimization problem
    std::vector<double> x{1,2};
    std::vector<double> g{1,2};
    double res; 
    res  = objective_coverage(x, g, nullptr);
    //assert((res -1.0) < eps);
    //assert((g[0] -1.0) < eps);

    // test constraint function
    const std::size_t n = 2;
    constraint_coeffs<n>  c;
    // xstar (2,2)
    c.coeffs << 2., 2.;
    c.sign = -1;
    // t=2, x1=1, x2=1
    std::vector<double> x_con{2,1,1};
    std::vector<double> g_con{1,1,1};
    double constraint_res;
    constraint_res = coverage_constraint<n>(x_con, g_con, &c); 
    //assert(abs(constraint_res) < eps);
    //-1* (((1,1) - (2,2))**2)  + 2 = 0
    c.sign = 1.;
    constraint_res = coverage_constraint<n>(x_con, g_con, &c); 
    //1* (((1,1) - (2,2))**2)  + 2 = 4 
   // assert(abs(constraint_res-4) < eps );

    // test optimization
    std::vector<double> ub{HUGE_VAL,4,3};
    std::vector<double> lb{10,0,0};
    Vector2d x_star_1;
    x_star_1 << 1,1;
    Vector2d x_star_2;
    x_star_2 << 2,2;
    opt cover_opt = opt("LD_MMA",n+1);
    opt cover_opt_local = opt("LD_SLSQP", n+1);
    cover_opt_local.set_xtol_rel(1e-4);
   // cover_opt.set_local_optimizer(cover_opt_local);
    cover_opt.set_upper_bounds(ub);
    cover_opt.set_lower_bounds(lb);


    constraint_coeffs<n> c1;
    c1.sign = 1;
    c1.coeffs = x_star_1;
    constraint_coeffs<n> c2;
    c2.coeffs = x_star_1;
    c2.sign = -1;
    constraint_coeffs<n> c3;
    c3.sign = 1;
    c3.coeffs = x_star_2;
    constraint_coeffs<n> c4;
    c4.coeffs = x_star_2;
    c4.sign = -1;
    cover_opt.add_inequality_constraint(coverage_constraint<n>, &c1, 1e-8);
    cover_opt.add_inequality_constraint(coverage_constraint<n>, &c1, 1e-8);
    cover_opt.add_inequality_constraint(coverage_constraint<n>, &c3, 1e-8);
    cover_opt.add_inequality_constraint(coverage_constraint<n>, &c4, 1e-8);
    cover_opt.set_xtol_rel(1e-4);
    cover_opt.set_min_objective(objective_coverage, NULL);
    cover_opt.set_maxeval(2000);
    double minf;
    std::vector<double> x0{20,4,0};
    cover_opt.optimize(x0, minf);
    std::cout << x0[0] << " " <<x0[1] <<" " << x0[2] << std::endl;
    std::cout << c1.coeffs <<std::endl;
}
    
