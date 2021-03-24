#include <iostream>
#include "optimizer.hpp"
#include <vector>
#include <nlopt.hpp>

using namespace nlopt;

int main()
{   

    // biased optimizer
    // default constructur
    
    const std::size_t n = 2;
    const int n_iter = 5;
    std::vector<double> lb{0,0};
    std::vector<double> ub{2,4};

    biased_optimizer<n> bopt; 
    // set bounds
    bopt.set_bounds(lb, ub);
    // add constraints
    constraint_coeffs<n> c;
    c.coeffs << 1., 1.;
    c.cons = 2.;
    c.type = "ineq";
    bopt.add_constraints(c);
    bopt.run(n_iter);
    // equality and inequality constraints
    biased_optimizer<n> bopt2; 
    // set bounds
    bopt2.set_bounds(lb, ub);
    // add constraints
    constraint_coeffs<n> eqc;
    eqc.coeffs << 1. , 1.;
    eqc.cons = 2.;
    eqc.type = "eq";
    constraint_coeffs<n> ineqc;
    ineqc.coeffs << 1. , 1.;
    ineqc.cons = 2.;
    ineqc.type = "ineq";
    bopt2.add_constraints(eqc, ineqc);
    bopt2.run(n_iter);
    // initialize with constraints and bounds
    biased_optimizer<n> bopt3(eqc, ineqc, lb, ub);
    bopt3.run(n_iter);
    // test run
    biased_optimizer<n> bopt4(ineqc, lb, ub);
    bopt4.run(n_iter);
    for (int i=0; i<n_iter; i++)
    {   
        for (std::vector<double>::const_iterator j = bopt4.results[i].begin(); j != bopt4.results[i].end(); ++j)
        {
            std::cout << *j << " ";
        }
        std::cout << "\n";
    }

    return 0;
}
