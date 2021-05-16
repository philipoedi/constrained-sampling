#include <iostream>
#include "optimizer.hpp"
#include <vector>
#include <nlopt.hpp>

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
    c.constype = "linear";
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
    eqc.constype = "linear";
    constraint_coeffs<n> ineqc;
    ineqc.coeffs << 1. , 1.;
    ineqc.cons = 2.;
    ineqc.type = "ineq";
    ineqc.constype = "linear";
    bopt2.add_constraints(eqc, ineqc);
    bopt2.run(n_iter);
    // initialize with constraints and bounds
    biased_optimizer<n> bopt3(eqc, ineqc, lb, ub);
    bopt3.run(n_iter);
    // test run
    biased_optimizer<n> bopt4(ineqc, lb, ub);
    bopt4.run(n_iter);
    std::vector<std::vector<double>> results;
    bopt4.results(results);
    print_matrix(results);  
    check_within_bounds(lb,ub,results);
    // test quadratic constraint
    constraint_coeffs<n> quad;
    double res;
    quad.P = Matrix2d::Identity();
    quad.type = "ineq";
    quad.constype = "quadratic";
    // test P term - identity
    std::vector<double> gradient(n);
    res = quadratic_constraint<n>(ub, gradient, &quad);
    std::cout << res <<std::endl;
    assert (res == 10.);
    quad.q(0) = 1;
    quad.q(1) = 1;  
    quad.r = 1;
    res = quadratic_constraint<n>(ub, gradient, &quad);
    std::cout << res <<std::endl;
    assert (res == 15);
    std::cout << gradient[0] << " " << gradient[1] << std::endl;
    Vector2d grad_true{3,5};
    // check gradient
    for (int i=0; i<n; i++)
    {
        assert (grad_true(i) == gradient[i]);
    }

    // biased optimizer with quadratic constraint
    biased_optimizer<n> bopt5(quad, lb, ub);
    bopt5.results(results);
    print_matrix(results); 
    check_within_bounds(lb,ub,results);

    // slacked optimizer

    const std::size_t m{1};
    const std::size_t l{0};
    
    slack_optimizer<n,m,l> sopt;
    sopt.set_bounds(lb, ub);
    sopt.add_constraints(ineqc);
    std::vector<double> x0{2,2};
    double slack;
    // find slack for inequality only
    slack = sopt.find_slack(x0, ineqc);
    std::cout << "slack: " << slack << std::endl;
    assert (slack == 2);


    slack_data<n,m,l> sl_data;
    init_slack_data(sl_data);
    // test initialization of slack data 
    // sl_data.a holds vector a.T @ x for objective function
    std::cout << "slack_data init: \n" << sl_data.a << "\nslack_data end" << std::endl;
    const std::size_t l2{1};
    slack_data<n,m,l2> sl_data2;
    init_slack_data(sl_data2);
    std::cout << "slack_data init: \n" << sl_data2.a << "\nslack_data end" << std::endl;
    double res2;
    res2 = slack_objective<n,m,l>(x_slacked, g_slacked, & sl_data);
    std::cout << res2 << std::endl;
    // constructor with ineq
    slack_optimizer<n,m,l> sopt2(ineqc, lb, ub); 
    sopt2.set_bounds(lb,ub);
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
