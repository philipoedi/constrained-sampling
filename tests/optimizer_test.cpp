#include <iostream>
#include "optimizer.hpp"
#include <vector>
#include <nlopt.hpp>
#include "constraints.hpp"

using namespace nlopt;


void print_matrix(std::vector<std::vector<double>> mat)
{
    for (int i=0; i<mat.size(); i++)
    {   
        for (std::vector<double>::const_iterator j = mat[i].begin(); j != mat[i].end(); ++j)
        {
            std::cout << *j << " ";
        }
        std::cout << "\n";
    }
}

void check_within_bounds(std::vector<double> lb, std::vector<double> ub, std::vector<std::vector<double>> mat)
{
    for (int i=0; i<mat.size(); i++)
    {   
        for (int j=0; j<mat[i].size(); ++j)
        {
            assert (mat[i][j] <= ub[j]);
            assert (mat[i][j] >= lb[j]);
        }
    }
}

int main()
{   

    // Biased optimizer
    // default constructur
    
    const std::size_t n = 2;
    const int n_iter = 5;
    std::vector<double> lb{0,0};
    std::vector<double> ub{2,4};

    BiasedOptimizer<n> bopt; 
    // set bounds
    bopt.setBounds(lb, ub);

    // add constraints
    ConstraintCoeffs<n> c;
    c.coeffs << 1., 1.;
    c.cons = 2.;
    c.type = "ineq";
    c.constype = "linear";
    bopt.addConstraints(c);
    bopt.run(n_iter);
    
    // equality and inequality constraints
    BiasedOptimizer<n> bopt2; 
    // set bounds
    bopt2.setBounds(lb, ub);
    // add constraints
    ConstraintCoeffs<n> eqc;
    eqc.coeffs << 1. , 1.;
    eqc.cons = 2.;
    eqc.type = "eq";
    eqc.constype = "linear";
    ConstraintCoeffs<n> ineqc;
    ineqc.coeffs << 1. , 1.;
    ineqc.cons = 2.;
    ineqc.type = "ineq";
    ineqc.constype = "linear";
    bopt2.addConstraints(eqc, ineqc);
    bopt2.run(n_iter);
    // initialize with constraints and bounds
    BiasedOptimizer<n> bopt3(eqc, ineqc, lb, ub);
    bopt3.run(n_iter);
    // test run
    BiasedOptimizer<n> bopt4(ineqc, lb, ub);
    bopt4.run(n_iter);
    std::vector<std::vector<double>> results;
    bopt4.results(results);
    print_matrix(results);  
    check_within_bounds(lb,ub,results);
    // test quadratic constraint
    ConstraintCoeffs<n> quad;
    double res;
    quad.P = Matrix2d::Identity();
    quad.type = "ineq";
    quad.constype = "quadratic";
    // test P term - identity
    std::vector<double> gradient(n);
    res = quadraticConstraint<n>(ub, gradient, &quad);
    std::cout << res <<std::endl;
    assert (res == 10.);
    quad.q(0) = 1;
    quad.q(1) = 1;  
    quad.r = 1;
    res = quadraticConstraint<n>(ub, gradient, &quad);
    std::cout << res <<std::endl;
    assert (res == 15);
    std::cout << gradient[0] << " " << gradient[1] << std::endl;
    Vector2d grad_true{3,5};
    // check gradient
    for (int i=0; i<n; i++)
    {
        assert (grad_true(i) == gradient[i]);
    }

    // Biased optimizer with quadratic constraint
    BiasedOptimizer<n> bopt5(quad, lb, ub);
    bopt5.results(results);
    print_matrix(results); 
    check_within_bounds(lb,ub,results);

    // slacked optimizer

    const std::size_t m{1};
    const std::size_t l{0};
    
    SlackOptimizer<n,m,l> sopt;
    sopt.setBounds(lb, ub);
    sopt.addConstraints(ineqc);
    std::vector<double> x0{2,2};
    double slack;
    // find slack for inequality only
    slack = sopt.findSlack(x0, ineqc);
    std::cout << "slack: " << slack << std::endl;
    assert (slack == 2);


    SlackData<n,m,l> sl_data;
    initSlackData(sl_data);
    // test initialization of slack data 
    // sl_data.a holds vector a.T @ x for objective function
    std::cout << "SlackData init: \n" << sl_data.a << "\nSlackData end" << std::endl;
    const std::size_t l2{1};
    SlackData<n,m,l2> sl_data2;
    initSlackData(sl_data2);
    std::cout << "SlackData init: \n" << sl_data2.a << "\nSlackData end" << std::endl;
    double res2;
    std::vector<double> x_slacked, g_slacked;
    res2 = slackObjective<n,m,l>(x_slacked, g_slacked, & sl_data);
    std::cout << res2 << std::endl;
    // constructor with ineq
    SlackOptimizer<n,m,l> sopt2(ineqc, lb, ub); 
    sopt2.setBounds(lb,ub);
    sopt2.sample(x_slacked);


    // check run
    sopt2.run(n_iter);
    std::vector<std::vector<double>> results_slack;
    std::cout << "checking results of slack optimizer" << std::endl;
    sopt2.results(results_slack);
    print_matrix(results_slack);
    std::vector<std::vector<double>> samples;
    std::cout << "cehcking samples of slack optimizer" << std::endl;
    sopt2.samples(samples);
    print_matrix(samples);


    return 0;


}
