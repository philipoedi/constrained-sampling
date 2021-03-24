#ifndef OPTIMIZER_H
#define OPTIMIZER_H
#include <iomanip>
#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <nlopt.hpp>
#include <list>
#include <string>
#include "sampler.hpp"
#include "utils.hpp"

using namespace nlopt;
using namespace Eigen;

template<std::size_t n>
struct bias {
    Matrix<double, n, 1> x0;
};

template<std::size_t n>
struct constraint_coeffs {
    Matrix<double, n, 1> coeffs = Matrix<double, n, 1>::Zero();
    Matrix<double, n, 1> q = Matrix<double, n, 1>::Zero();
    Matrix<double, n, n> P = Matrix<double, n, n>::Zero();
    double cons;
    double sign;
    double r{0};
    std::string type;
};

template<std::size_t n> double biased_objective(const std::vector<double>& x, std::vector<double>& grad, void*data )
{
    typedef Matrix<double, n, 1> vec;
    bias<n> *b = (bias<n>*) data;
    vec x_vec(x.data());
    vec x_x0 = x_vec - b->x0;
    if(!grad.empty()){       
        for (std::size_t i = 0; i<x_x0.size() ;++i){
            grad[i] = 2*x_x0[i];
        }
    }   
    return x_x0.transpose()*x_x0;    
}

template<std::size_t n> double linear_constraint(const std::vector<double>& x, std::vector<double> &grad, void*data)
{
    typedef Matrix<double, n, 1> vec;
    constraint_coeffs<n> *c = (constraint_coeffs<n>*) data;
    vec x_vec(x.data());
    if (!grad.empty()){
        utils::copy_eig2vec(c->coeffs, grad);
    }
    return x_vec.transpose() * c->coeffs - c->cons;
}


template<std::size_t n> double quadratic_constraint(const std::vector<double> x, std::vector<double>& grad, void* data)
{
    typedef Matrix<double, n, 1> vec;
    double res = 0;
    constraint_coeffs<n> *c = (constraint_coeffs<n>*) data;
    vec x_vec(x.data());
    // 0.5 * x.T@P@x + q.T@x+ r
    if (!grad.empty()){
        utils::copy_eig2vec(c->P.transpose()*x_vec + c->q, grad);
    }
    res += 0.5 * x_vec.transpose() * c->P * x_vec; 
    res += x_vec.transpose()*c->q;  
    return res + c->r;
}


template<std::size_t n>
class biased_optimizer
{
    typedef Matrix<double,n,1> Vector;

    public:

        biased_optimizer();
        biased_optimizer(
            constraint_coeffs<n>& eqcons,
            constraint_coeffs<n>& ineqcons, 
            const std::vector<double>& lb, 
            const std::vector<double>& ub); 
        biased_optimizer(
            constraint_coeffs<n>& cons,
            const std::vector<double>& lb, 
            const std::vector<double>& ub); 
        void run(const int niter);
        void save();
        void add_constraints(
            constraint_coeffs<n>& eqcons,
            constraint_coeffs<n>& ineqcons);
        void add_constraints(constraint_coeffs<n>& cons);
        void set_bounds(const std::vector<double>& lb, const std::vector<double>& ub);

        std::vector<std::vector<double>> results;
        std::vector<std::vector<double>> samples;
    private:
        
        opt opt_{"AUGLAG_EQ",n};
        opt local_opt_{"LD_SLSQP",n};
        double minf_;
        uniform_sampler<n> uni_;
        bias<n> b_;
};


template<std::size_t n>
biased_optimizer<n>::biased_optimizer()
{
    local_opt_.set_xtol_rel(1e-4);
    opt_.set_xtol_rel(1e-4);
    opt_.set_local_optimizer(local_opt_);
}

template<std::size_t n>
biased_optimizer<n>::biased_optimizer(
    constraint_coeffs<n>& eqcons,
    constraint_coeffs<n>& ineqcons,
    const std::vector<double>& lb,
    const std::vector<double>& ub)
{
    local_opt_.set_xtol_rel(1e-4);
    opt_.set_xtol_rel(1e-4);
    opt_.set_local_optimizer(local_opt_);
    this->set_bounds(lb, ub);
    this->add_constraints(eqcons, ineqcons);
}

template<std::size_t n>
biased_optimizer<n>::biased_optimizer(
    constraint_coeffs<n>& cons,
    const std::vector<double>& lb,
    const std::vector<double>& ub)
{
    local_opt_.set_xtol_rel(1e-4);
    opt_.set_xtol_rel(1e-4);
    opt_.set_local_optimizer(local_opt_);
    this->set_bounds(lb, ub);
    this->add_constraints(cons);
}

template<std::size_t n>
void biased_optimizer<n>::set_bounds(const std::vector<double>& lb, const std::vector<double>& ub)
{
    uni_.set_bounds(lb,ub);
    opt_.set_upper_bounds(ub);
    opt_.set_lower_bounds(lb);
}

template<std::size_t n>
void biased_optimizer<n>::add_constraints(
    constraint_coeffs<n>& eqcons,
    constraint_coeffs<n>& ineqcons)
{
    opt_.add_equality_constraint(linear_constraint<n>, &eqcons, 1e-8);
    opt_.add_inequality_constraint(linear_constraint<n>, &ineqcons, 1e-8);
}

template<std::size_t n>
void biased_optimizer<n>::add_constraints(constraint_coeffs<n>& cons)
{
    if (cons.type == "eq")
    {
        opt_.add_equality_constraint(linear_constraint<n>, &cons, 1e-8);
    }
    else if (cons.type == "ineq")
    {
        opt_.add_inequality_constraint(linear_constraint<n>, &cons, 1e-8);
    }
}


template<std::size_t n>
void biased_optimizer<n>::run(const int niter)
{
    std::vector<double> x(n);
    results.resize(niter);
    samples.resize(niter);
    for (int i=0; i<niter; ++i)
    {
        b_.x0 = uni_.sample();
        utils::copy_eig2vec(b_.x0, x);
        opt_.set_min_objective(biased_objective<n>, &b_);
        try
        {
            opt_.optimize(x, minf_);
            results[i].resize(n);
            samples[i].resize(n);
            utils::copy_vec2vec(x, results[i]);
            utils::copy_eig2vec(b_.x0, samples[i]);
        }
        catch(std::exception &e)
        {
            std::cerr << e.what() << std::endl;
        }
        
    }
}


#endif



