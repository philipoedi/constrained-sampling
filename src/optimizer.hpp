#include <iomanip>
#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <nlopt.hpp>

using namespace nlopt;
using namespace Eigen;

template<std::size_t n>
struct bias {
    Matrix<double, n, 1> x0;
};

template<std::size_t n>
struct constraint_coeffs {
    Matrix<double, n, 1> coeffs;
    double cons;
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
    vec coeffs(c->coeffs.data());
    if (!grad.empty()){
        for (std::size_t i = 0; i<x_vec.size(); ++i){
           grad[i] = coeffs[i];
        }
    }
    return x_vec.transpose() * c->coeffs - c->cons;
}

